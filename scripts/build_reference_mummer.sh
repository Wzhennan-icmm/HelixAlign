#!/usr/bin/env bash
# Build the unmodified, checksum-pinned official exact-matching executable.
set -euo pipefail
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 BUILD_DIRECTORY" >&2
  exit 2
fi
mkdir -p "$1"
MUMMER_BUILD_ROOT=$(cd "$1" && pwd)
MUMMER_ARCHIVE="$MUMMER_BUILD_ROOT/mummer-4.0.1.tar.gz"
MUMMER_SOURCE="$MUMMER_BUILD_ROOT/mummer-4.0.1"
if [[ -e "$MUMMER_SOURCE" ]]; then
  echo "Refusing to reuse an existing source tree: $MUMMER_SOURCE" >&2
  exit 2
fi
curl -L --fail --retry 2 \
  https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz \
  -o "$MUMMER_ARCHIVE"
printf '%s  %s\n' bc20ae2701a0b2e323e4e515b7cfa18a0f0cb34a4ff5844b289b2de0154e3d3e "$MUMMER_ARCHIVE" | sha256sum -c -
tar --no-same-owner -xzf "$MUMMER_ARCHIVE" -C "$MUMMER_BUILD_ROOT"
cd "$MUMMER_SOURCE"
./configure --prefix="$MUMMER_BUILD_ROOT/install"
make -j "${MUMMER_BUILD_JOBS:-4}" mummer
./mummer -version
printf 'Official executable (retain build tree): %s/mummer\n' "$MUMMER_SOURCE"
