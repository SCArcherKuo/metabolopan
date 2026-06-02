use std::time::SystemTime;

use metabolopan::logging::{LogLayer, LogLine, LogStore};
use tracing::{Level, info};
use tracing_subscriber::Registry;
use tracing_subscriber::layer::SubscriberExt;

#[test]
fn bounded_buffer_drops_oldest() {
    let store = LogStore::new(3);
    for i in 0..5 {
        store.push(LogLine {
            level: Level::INFO,
            timestamp: SystemTime::now(),
            target: "test".to_string(),
            message: format!("event {i}"),
        });
    }
    let snap = store.snapshot();
    assert_eq!(snap.len(), 3);
    let messages: Vec<&str> = snap.iter().map(|l| l.message.as_str()).collect();
    assert_eq!(messages, vec!["event 2", "event 3", "event 4"]);
}

#[test]
fn layer_routes_events_to_store() {
    let store = LogStore::new(100);
    let subscriber = Registry::default().with(LogLayer::new(store.clone()));

    tracing::subscriber::with_default(subscriber, || {
        info!("hello from layer test");
        info!(key = "value", "structured event");
    });

    let snap = store.snapshot();
    assert!(
        snap.iter()
            .any(|l| l.message.contains("hello from layer test")),
        "expected first message in store; got: {:?}",
        snap.iter().map(|l| &l.message).collect::<Vec<_>>()
    );
    assert!(
        snap.iter().any(|l| l.message.contains("structured event")),
        "expected second message in store"
    );
}

#[test]
fn clear_empties_buffer() {
    let store = LogStore::new(10);
    store.push(LogLine {
        level: Level::INFO,
        timestamp: SystemTime::now(),
        target: "test".to_string(),
        message: "x".to_string(),
    });
    assert!(!store.is_empty());
    store.clear();
    assert!(store.is_empty());
}
