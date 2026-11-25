#!/bin/bash
# T1-B FASE 3: bounded post-processing (READ-ONLY on archives; 1 core, 15 min).
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
ulimit -v 2457600
cd "$DIR"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
timeout -s TERM --kill-after=5 900 /usr/bin/time -v taskset -c 0 \
  python3 -I t1b_validate.py > validate.log 2> validate_resource_report.txt
rc=$?
echo "exit_code=$rc" >> validate_resource_report.txt
exit $rc
