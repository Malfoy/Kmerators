#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cargo build \
  --manifest-path "$ROOT/Cargo.toml" \
  --release \
  --target wasm32-unknown-unknown

mkdir -p "$ROOT/web/wasm"
cp \
  "$ROOT/target/wasm32-unknown-unknown/release/kmerators_wasm_core.wasm" \
  "$ROOT/web/wasm/kmerators_wasm_core.wasm"

ls -lh "$ROOT/web/wasm/kmerators_wasm_core.wasm"
