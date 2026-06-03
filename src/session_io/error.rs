//! Error types for the `session-settings-io` capability.
//!
//! Five variants distinguish the failure modes a save / load / hash-diff
//! operation can hit. The `Display` impls produce single-line user-readable
//! strings suitable for the error-toast modal in `app-shell`; the OS-level
//! cause stays in the `source` chain for developers.
//!
//! See the `session-settings-io` capability spec (Requirement:
//! `SnapshotError` SHALL distinguish IO, parse, version, and write
//! failures) for the normative contract.

use std::fmt;
use std::io;
use std::path::PathBuf;

#[derive(Debug)]
pub enum SnapshotError {
    /// File read failed during load (file not found, permission denied,
    /// IO error before any JSON parsing happened).
    Io { path: PathBuf, source: io::Error },

    /// `serde_json` parse failure during load. `line` / `column` come from
    /// `serde_json::Error::line()` / `column()`. Either field is `0` when
    /// serde_json could not locate the error (typically an unexpected EOF
    /// on binary input); in that case the `Display` rendering omits the
    /// parenthetical position to avoid the misleading `(line 0 column 0)`.
    JsonParse {
        path: PathBuf,
        line: usize,
        column: usize,
        message: String,
    },

    /// The file's `schema_version` field is not equal to `SCHEMA_VERSION`.
    /// Checked BEFORE the body deserialise so a future-version payload
    /// cannot accidentally parse under the current contract.
    UnsupportedSchemaVersion { found: u32, expected: u32 },

    /// SHA-256 computation failed during save or `diff_input_hashes`.
    /// Distinct from `Io` so the caller can surface "we couldn't verify
    /// inputs" separately from "we couldn't read the snapshot file".
    HashIo { path: PathBuf, source: io::Error },

    /// JSON serialise OR file write failed during save. The function MUST
    /// NOT leave a partial file at `path` when returning this variant.
    WriteIo { path: PathBuf, source: io::Error },
}

impl fmt::Display for SnapshotError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SnapshotError::Io { path, source } => {
                write!(f, "Cannot read file {}: {}", path.display(), source)
            }
            SnapshotError::JsonParse {
                path,
                line,
                column,
                message,
            } => {
                if *line == 0 && *column == 0 {
                    write!(
                        f,
                        "Settings file {} is not valid JSON: {}",
                        path.display(),
                        message
                    )
                } else {
                    write!(
                        f,
                        "Settings file {} is not valid JSON (line {} column {}): {}",
                        path.display(),
                        line,
                        column,
                        message
                    )
                }
            }
            SnapshotError::UnsupportedSchemaVersion { found, expected } => {
                write!(
                    f,
                    "This settings file uses schema version {}; this app expects version {}. The file was likely saved by a different app version.",
                    found, expected
                )
            }
            SnapshotError::HashIo { path, source } => {
                write!(
                    f,
                    "Cannot verify input file hash for {}: {}",
                    path.display(),
                    source
                )
            }
            SnapshotError::WriteIo { path, source } => {
                write!(f, "Cannot save settings to {}: {}", path.display(), source)
            }
        }
    }
}

impl std::error::Error for SnapshotError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SnapshotError::Io { source, .. }
            | SnapshotError::HashIo { source, .. }
            | SnapshotError::WriteIo { source, .. } => Some(source),
            SnapshotError::JsonParse { .. } | SnapshotError::UnsupportedSchemaVersion { .. } => {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn io_err() -> io::Error {
        io::Error::new(io::ErrorKind::NotFound, "no such file")
    }

    #[test]
    fn display_io() {
        let e = SnapshotError::Io {
            path: PathBuf::from("/tmp/x.json"),
            source: io_err(),
        };
        assert_eq!(e.to_string(), "Cannot read file /tmp/x.json: no such file");
    }

    #[test]
    fn display_json_parse_with_location() {
        let e = SnapshotError::JsonParse {
            path: PathBuf::from("/tmp/x.json"),
            line: 7,
            column: 15,
            message: "missing field `settings`".to_string(),
        };
        assert_eq!(
            e.to_string(),
            "Settings file /tmp/x.json is not valid JSON (line 7 column 15): missing field `settings`"
        );
    }

    #[test]
    fn display_json_parse_without_location() {
        let e = SnapshotError::JsonParse {
            path: PathBuf::from("/tmp/x.json"),
            line: 0,
            column: 0,
            message: "EOF while parsing".to_string(),
        };
        assert_eq!(
            e.to_string(),
            "Settings file /tmp/x.json is not valid JSON: EOF while parsing"
        );
    }

    #[test]
    fn display_unsupported_schema_version() {
        let e = SnapshotError::UnsupportedSchemaVersion {
            found: 2,
            expected: 1,
        };
        assert_eq!(
            e.to_string(),
            "This settings file uses schema version 2; this app expects version 1. The file was likely saved by a different app version."
        );
    }

    #[test]
    fn display_hash_io() {
        let e = SnapshotError::HashIo {
            path: PathBuf::from("/data/POS.txt"),
            source: io_err(),
        };
        assert_eq!(
            e.to_string(),
            "Cannot verify input file hash for /data/POS.txt: no such file"
        );
    }

    #[test]
    fn display_write_io() {
        let e = SnapshotError::WriteIo {
            path: PathBuf::from("/tmp/x.json"),
            source: io_err(),
        };
        assert_eq!(
            e.to_string(),
            "Cannot save settings to /tmp/x.json: no such file"
        );
    }
}
