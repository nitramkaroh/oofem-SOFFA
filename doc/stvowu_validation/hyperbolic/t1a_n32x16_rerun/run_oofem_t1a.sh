#!/bin/bash
# T1-A: single OOFEM re-run of the archived n32x16 case (authorized 2026-10-27).
# Envelope: 1 process, 1 core (taskset -c 0), virtual-memory cap 8 GB (ulimit -v),
# hard timeout 60 s (timeout -s TERM --kill-after=5).
# Invocation identical to the original (run_coupled_independent.py):
#   oofem -f <folder>/case.in   with cwd = <folder>
set -u
DIR="$(cd "$(dirname "$0")/oofem" && pwd)"
OOFEM=/home/marco/Documents/Work/oofem-SOFFA/build/debug_magneto/oofem
ulimit -v 8388608
cd "$DIR"
timeout -s TERM --kill-after=5 60 /usr/bin/time -v taskset -c 0 "$OOFEM" -f "$DIR/case.in" > run.log 2> resource_report.txt
rc=$?
echo "exit_code=$rc" >> resource_report.txt
exit $rc
