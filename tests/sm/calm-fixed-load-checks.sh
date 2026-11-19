#
# CALM (solverType "calm") when the linearised arc-length constraint has no solution: calm-fixed-load.in
# with the arc length 0.1, ten times shorter than the displacement caused by the fixed load. The corrector
# line then misses the constraint circle (negative discriminant) after every restart, so the solver must
# restart the step maxrestarts times and then report a failed step to the engineering model, which stops
# here because no step reduction strategy is set. It must not stop inside the solver with an error.
# Run by ctest as: bash calm-fixed-load-checks.sh <oofem executable>
#
OOFEM=$1
case "$OOFEM" in
    /* | ?:*) ;;
    *) OOFEM="$PWD/$OOFEM" ;;
esac
INPUT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) || exit 1
INPUT="$INPUT_DIR/calm-fixed-load.in"
WORK=$(mktemp -d) || exit 1
trap 'cd / && rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1

awk '
    /^#%BEGIN_CHECK%/ { exit }
    /^errorcheck[[:space:]]*$/ { next }
    NR == 1 { print "no-intersection.out"; next }
    {
        gsub("stepLength 1.5", "stepLength 0.1")
        gsub("nmodules 1", "nmodules 0")
        print
    }
' "$INPUT" > no-intersection.in

if "$OOFEM" -f no-intersection.in > no-intersection.log 2>&1; then
    echo 'FAIL no-intersection: the analysis accepted a step without a solution of the constraint'
    cat no-intersection.log
    exit 1
fi
restarts=$(grep -c 'discriminant is negative' no-intersection.log)
if [ "$restarts" -eq 3 ] &&
   grep -q 'arc-length constraint could not be satisfied after 2 restarts' no-intersection.log &&
   grep -q 'Nonlinear solver failed to converge at step 1' no-intersection.log &&
   ! grep -q "can't continue further" no-intersection.log; then
    echo 'OK   no-intersection: 2 restarts, then the failed step reached the engineering model'
else
    echo "FAIL no-intersection: expected 3 attempts and a failed step, found $restarts attempts"
    cat no-intersection.log
    exit 1
fi
