#!/usr/bin/env bash
# Starts the local Hugo dev server with live reload.
# Usage: ./start.sh [extra hugo server flags]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# In case Hugo was installed to ~/.local/bin (no sudo) rather than a
# system-wide location already on PATH.
export PATH="$HOME/.local/bin:$PATH"

if ! command -v hugo >/dev/null 2>&1; then
  echo "error: hugo not found on PATH." >&2
  echo "Install it from https://github.com/gohugoio/hugo/releases" >&2
  echo "(grab hugo_extended_<version>_<os>-<arch>.tar.gz, extract the" >&2
  echo "'hugo' binary, and put it somewhere on your PATH, e.g. ~/.local/bin)." >&2
  exit 1
fi

echo "Starting Hugo dev server — visit http://localhost:1313"
exec hugo server -D "$@"
