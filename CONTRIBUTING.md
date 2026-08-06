# Contributing to Metabolopan

Thanks for your interest in contributing! Metabolopan is a single-binary Rust
desktop app (`eframe`/`egui` GUI + `tokio` for network I/O) that runs the whole
MS-DIAL → DAM → KEGG ORA pipeline in-process. This guide covers dev setup and
code/commit conventions.

## Dev setup

- **Rust 1.85+** (the crate is on the Rust 2024 edition).
- Run the GUI: `cargo run --release` (debug builds are noticeably slow on a
  ~13k-feature DAM run).
- Type-check: `cargo check`
- Lint: `cargo clippy --all-targets -- -D warnings` (must exit 0)
- Format: `cargo fmt --all` (Rust 2024 style)
- Tests: `cargo test` (800+ tests; HTTP integration tests use `wiremock` and do
  not hit the network). The KEGG/HTTP smoke tests are `#[ignore]`d by default —
  run them explicitly with `cargo test -- --ignored` only when you intend to
  reach the live KEGG/PubChem services.

## Proposing changes

- **Feature-level work** — open an issue first to discuss the design and scope,
  then submit a PR referencing it.
- **Bug fixes / typos / small refactors** — go straight to a PR.

Explain the **why** in the PR description. Where a decision constrains the code
— an invariant, an ordering, a value that must not change — state that
constraint in a doc-comment too: a PR description is not in the tree a later
reader has.

## Code style

- `cargo fmt` (enforced) and `cargo clippy --all-targets -- -D warnings`
  (enforced) both pass before you push.
- Keep comments as short as the point allows, and no shorter: a comment that
  states a constraint the code must satisfy earns the lines it takes.
- **UI text is ASCII** — the default egui font has gaps in its Unicode coverage
  (e.g. it lacks `←` and `›`, which is why Back/stepper separators use `<` and
  `>`). Keep on-screen strings ASCII unless you've verified the glyph renders.
- Reference theme colors by name from `src/theme.rs` — no inline
  `Color32::from_rgb(...)` literals under `src/app.rs` or `src/ui/`.

## Test conventions

- Pure functions → in-module `#[cfg(test)] mod tests`.
- Cross-module / integration → `tests/<name>_test.rs`.
- HTTP integration → `wiremock` (never hit the real network in CI).
- End-to-end GUI tests are not yet automated; UI changes rely on manual
  verification — include before/after screenshots in the PR.

## Commit conventions

- Imperative subject line, ≤ 70 chars.
- Body explains **why**, not what.
- Don't amend or force-push commits that are already pushed; fix forward with a
  new commit.

## Reporting bugs

The GUI's log pane has a **[Download bug report…]** button that produces a
privacy-redacted zip (logs + app-state snapshot + input/cache summaries; never
your raw MS-DIAL/metadata files, never your home path). Attaching that zip to a
bug report is the fastest way to get a fix.

## License

By contributing, you agree that your contributions are licensed under the
project's [Apache License 2.0](LICENSE).
