use std::collections::VecDeque;
use std::fmt;
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

use tracing::Level;
use tracing::field::{Field, Visit};
use tracing_subscriber::Layer;
use tracing_subscriber::layer::Context;
use tracing_subscriber::registry::LookupSpan;

#[derive(Debug, Clone)]
pub struct LogLine {
    pub level: Level,
    pub timestamp: SystemTime,
    pub target: String,
    pub message: String,
}

#[derive(Debug, Clone)]
pub struct LogStore {
    inner: Arc<Mutex<VecDeque<LogLine>>>,
    cap: usize,
}

impl LogStore {
    pub fn new(cap: usize) -> Self {
        Self {
            inner: Arc::new(Mutex::new(VecDeque::with_capacity(cap.min(1024)))),
            cap,
        }
    }

    pub fn push(&self, line: LogLine) {
        let mut buf = self.inner.lock().expect("log store poisoned");
        if buf.len() == self.cap {
            buf.pop_front();
        }
        buf.push_back(line);
    }

    pub fn snapshot(&self) -> Vec<LogLine> {
        let buf = self.inner.lock().expect("log store poisoned");
        buf.iter().cloned().collect()
    }

    pub fn clear(&self) {
        let mut buf = self.inner.lock().expect("log store poisoned");
        buf.clear();
    }

    pub fn len(&self) -> usize {
        self.inner.lock().expect("log store poisoned").len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct LogLayer {
    store: LogStore,
}

impl LogLayer {
    pub fn new(store: LogStore) -> Self {
        Self { store }
    }
}

#[derive(Default)]
struct MessageVisitor {
    message: String,
}

impl Visit for MessageVisitor {
    fn record_str(&mut self, field: &Field, value: &str) {
        if field.name() == "message" {
            if !self.message.is_empty() {
                self.message.push(' ');
            }
            self.message.push_str(value);
        } else {
            if !self.message.is_empty() {
                self.message.push(' ');
            }
            self.message.push_str(&format!("{}={value}", field.name()));
        }
    }

    fn record_debug(&mut self, field: &Field, value: &dyn fmt::Debug) {
        if field.name() == "message" {
            if !self.message.is_empty() {
                self.message.push(' ');
            }
            self.message.push_str(&format!("{value:?}"));
        } else {
            if !self.message.is_empty() {
                self.message.push(' ');
            }
            self.message
                .push_str(&format!("{}={value:?}", field.name()));
        }
    }
}

impl<S> Layer<S> for LogLayer
where
    S: tracing::Subscriber + for<'a> LookupSpan<'a>,
{
    fn on_event(&self, event: &tracing::Event<'_>, _ctx: Context<'_, S>) {
        let metadata = event.metadata();
        let mut visitor = MessageVisitor::default();
        event.record(&mut visitor);
        self.store.push(LogLine {
            level: *metadata.level(),
            timestamp: SystemTime::now(),
            target: metadata.target().to_string(),
            message: visitor.message,
        });
    }
}
