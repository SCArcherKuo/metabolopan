use std::sync::Mutex;

use chrono::Utc;
use eframe::NativeOptions;
use egui::ViewportBuilder;
use tokio::runtime::Builder as TokioRuntimeBuilder;
use tracing_subscriber::{
    EnvFilter, Layer, Registry, fmt, layer::SubscriberExt, util::SubscriberInitExt,
};

use metabolopan::app::App;
use metabolopan::diagnostics::{
    CompactUtcTime, clear_stale_session_logs, file_sink_env_filter, session_log_dir,
    session_log_path, try_open_session_log,
};
use metabolopan::logging::{LogLayer, LogStore};

const STALE_LOG_MAX_AGE_DAYS: u32 = 7;

/// Decode the embedded placeholder window icon (`assets/icon.png`) into the
/// RGBA form egui wants for the runtime window / Dock / taskbar icon.
///
/// The PNG is `include_bytes!`-baked into the binary so there is no external
/// file dependency at run time (matches the single-binary architecture).
/// Replace `assets/icon.png` with real artwork to change the icon.
fn load_window_icon() -> egui::IconData {
    let image = image::load_from_memory(include_bytes!("../assets/icon.png"))
        .expect("embedded window icon PNG should decode")
        .to_rgba8();
    let (width, height) = (image.width(), image.height());
    egui::IconData {
        rgba: image.into_raw(),
        width,
        height,
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Resolve env filter directive (RUST_LOG, default "info").
    let directive = std::env::var("RUST_LOG").unwrap_or_else(|_| "info".to_string());
    let env_filter = EnvFilter::try_new(&directive).unwrap_or_else(|_| EnvFilter::new("info"));

    let log_store = LogStore::new(5000);

    // Per-session file sink setup. All failures here surface as deferred
    // WARN events emitted AFTER the subscriber is initialised (we have no
    // sink to emit to yet).
    let started_at = Utc::now();
    let pid = std::process::id();
    let log_dir = session_log_dir();
    let cleanup_report = clear_stale_session_logs(&log_dir, started_at, STALE_LOG_MAX_AGE_DAYS);
    let attempted_log_path = log_dir.join(session_log_path(started_at, pid));
    let (file_handle, session_log_path_resolved, file_sink_error) =
        match try_open_session_log(&log_dir, started_at, pid) {
            Ok((path, file)) => (Some(file), Some(path), None),
            Err(e) => (None, None, Some(e.to_string())),
        };

    let file_layer = file_handle.map(|f| {
        fmt::layer()
            .with_writer(Mutex::new(f))
            .with_ansi(false)
            .with_target(true)
            .with_timer(CompactUtcTime)
            .with_filter(file_sink_env_filter())
    });

    Registry::default()
        .with(env_filter)
        .with(fmt::layer().with_writer(std::io::stderr))
        .with(LogLayer::new(log_store.clone()))
        .with(file_layer)
        .init();

    // Deferred startup events (cleanup summary, file-sink failure WARN).
    tracing::info!(
        deleted = cleanup_report.deleted.len(),
        retained = cleanup_report.retained,
        skipped_missing_dir = cleanup_report.skipped_missing_dir,
        dir = %log_dir.display(),
        "session log cleanup completed"
    );
    for (path, err) in &cleanup_report.failures {
        tracing::warn!(
            path = %path.display(),
            error = %err,
            "failed to delete stale session log; continuing startup"
        );
    }
    if let Some(err) = file_sink_error {
        tracing::warn!(
            error = %err,
            path = %attempted_log_path.display(),
            "file sink unavailable; bug-report logs.txt will be a stub"
        );
    }

    tracing::info!(
        version = env!("CARGO_PKG_VERSION"),
        git_sha = option_env!("GIT_SHA").unwrap_or("unknown"),
        "metabolopan starting"
    );
    tracing::info!(
        cache_dir = %metabolopan::kegg::cache::cache_dir().display(),
        log_dir = %log_dir.display(),
        "runtime data directories"
    );

    let rt = TokioRuntimeBuilder::new_multi_thread()
        .enable_all()
        .thread_name("kegg-worker")
        .build()?;

    let native_options = NativeOptions {
        viewport: ViewportBuilder::default()
            .with_inner_size([1024.0, 720.0])
            .with_icon(load_window_icon()),
        ..Default::default()
    };

    eframe::run_native(
        "Metabolopan — Metabolomic Enrichment Analysis",
        native_options,
        Box::new(move |cc| {
            egui_extras::install_image_loaders(&cc.egui_ctx);
            metabolopan::theme::install(&cc.egui_ctx);
            Ok(Box::new(App::new(
                log_store,
                directive,
                rt,
                session_log_path_resolved,
            )))
        }),
    )?;
    Ok(())
}
