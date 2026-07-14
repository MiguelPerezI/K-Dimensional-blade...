#!/usr/bin/env bash
# Launch main15_gpu on the NVIDIA Tesla T4 OpenGL backend.
# Without __GLX_VENDOR_LIBRARY_NAME=nvidia the app falls back to Mesa llvmpipe
# (CPU software rasterization), which is what made the original main15 slow.
set -euo pipefail
cd "$(dirname "$0")"
export DISPLAY="${DISPLAY:-:0.0}"
export __GLX_VENDOR_LIBRARY_NAME=nvidia
exec ./main17_gpu "$@"
