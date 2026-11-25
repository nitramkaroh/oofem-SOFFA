#!/bin/bash
# T1-B FASE 2: single full DOP853 n128x64 integration — 1 process, 1 core,
# HARD 45-min outer timeout (2700 s), 2.4 GB virtual memory cap, in-process
# soft wall cap 40 min + NaN/Inf/growth guards inside t1b_dop853.py.
# NO second automatic attempt: this wrapper runs the process exactly once.
set -u
DIR="$(cd "$(dirname "$0")" && pwd)"
ulimit -v 2457600
cd "$DIR"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
timeout -s TERM --kill-after=5 2700 /usr/bin/time -v taskset -c 0 \
  python3 -I t1b_dop853.py full > run.log 2> resource_report.txt
rc=$?
echo "exit_code=$rc" >> resource_report.txt
exit $rc
