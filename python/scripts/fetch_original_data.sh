#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"
if [ -d data ]; then
  echo "data/ already exists; leaving it unchanged."
  exit 0
fi
git clone --depth 1 https://github.com/Reza219/Incremental_SLAM /tmp/Incremental_SLAM_py_port_upstream
cp -r /tmp/Incremental_SLAM_py_port_upstream/data ./data
echo "Copied upstream data/ into $(pwd)/data"
