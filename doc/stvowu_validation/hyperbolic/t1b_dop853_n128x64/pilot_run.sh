#!/bin/bash
# T1-B FASE 1: bounded pilot — 1 process, 1 core, HARD 60 s outer timeout.
# Memory cap: 2.4 GB virtual (baseline ~0.3 GB + 2 GB budget; measured need < 0.5 GB).
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
ulimit -v 2457600
cd "$DIR"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
timeout -s TERM --kill-after=5 60 /usr/bin/time -v taskset -c 0 \
  python3 -I t1b_dop853.py pilot > pilot.log 2> pilot_resource_report.txt
rc=$?
echo "exit_code=$rc" >> pilot_resource_report.txt
exit $rc
