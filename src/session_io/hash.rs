//! SHA-256 streaming helper for the `session-settings-io` capability.
//!
//! See `openspec/specs/session-settings-io/spec.md` (Requirement: SHA-256
//! helper SHALL stream files in chunks).

use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use sha2::{Digest, Sha256};

const CHUNK_BYTES: usize = 64 * 1024;

/// Compute the SHA-256 of a file's contents and return the digest as a
/// 64-character lowercase hex string. Streams the file in 64 KiB chunks
/// so multi-MB inputs do not need to be slurped into memory.
pub fn sha256_file(path: &Path) -> io::Result<String> {
    let mut file = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = [0u8; CHUNK_BYTES];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    let digest = hasher.finalize();
    Ok(hex_lower(&digest))
}

fn hex_lower(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len() * 2);
    for b in bytes {
        out.push_str(&format!("{:02x}", b));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_temp(contents: &[u8]) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("create tempfile");
        f.write_all(contents).expect("write tempfile");
        f.flush().expect("flush tempfile");
        f
    }

    #[test]
    fn sha256_abc_matches_rfc_vector() {
        let f = write_temp(b"abc");
        let digest = sha256_file(f.path()).expect("hash abc");
        assert_eq!(
            digest,
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn sha256_empty_file_matches_canonical_empty_digest() {
        let f = write_temp(b"");
        let digest = sha256_file(f.path()).expect("hash empty");
        assert_eq!(
            digest,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn sha256_two_mib_file_hashes_without_panic() {
        // 2 MiB of arbitrary bytes — exercises the chunked-read loop.
        let bytes: Vec<u8> = (0..(2 * 1024 * 1024)).map(|i| (i % 256) as u8).collect();
        let f = write_temp(&bytes);
        let digest = sha256_file(f.path()).expect("hash 2MiB");
        // Determinism check: hashing twice yields the same string.
        let digest2 = sha256_file(f.path()).expect("hash 2MiB again");
        assert_eq!(digest.len(), 64);
        assert_eq!(digest, digest2);
    }

    #[test]
    fn sha256_missing_file_returns_io_err() {
        let result = sha256_file(Path::new("/nonexistent/path/that/should/not/exist.bin"));
        assert!(result.is_err());
    }
}
