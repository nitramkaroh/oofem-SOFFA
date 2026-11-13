#
# Input checks of the NURBS and T-spline weights: one weight per control point, finite and positive
# (the rational basis functions are divided by sum(N*w)). The inputs are the NURBS patch test
# iga-nurbs-weights.in and the T-spline test iga-tspline-nip.in with modified weights; they are written
# into a temporary directory, so nothing is written to tests/sm.
# Run by ctest as: bash iga-weight-checks.sh <oofem executable>
#
OOFEM=$1
case "$OOFEM" in
    /* | ?:*) ;;
    *) OOFEM="$PWD/$OOFEM" ;;
esac
echo "target executable: $OOFEM"
SRC=$(cd "$(dirname "$0")" && pwd)

WORK=$(mktemp -d) || exit 1
trap 'cd / && rm -rf "$WORK"' EXIT
cd "$WORK" || exit 1

failures=0

# $1 case name, $2 source input, $3 its weights, $4 the modified weights, $5 expected part of the error
expect_error() {
    sed "s/$3/$4/" "$SRC/$2" > "$1.in"
    if "$OOFEM" -f "$1.in" > "$1.log" 2>&1; then
        echo "FAIL $1: must be rejected, but the analysis finished"
        failures=$((failures + 1))
    elif grep -q "$5" "$1.log"; then
        echo "OK   $1: rejected"
    else
        echo "FAIL $1: failed for another reason"
        cat "$1.log"
        failures=$((failures + 1))
    fi
}

N='weights 9 1.0 0.7 1.0 0.8 1.3 0.9 1.0 0.6 1.2'
T='weights 16 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'
COUNT='number of weights must equal the number of control points'
POSITIVE='weights must be finite and positive'

expect_error nurbs-count      iga-nurbs-weights.in "$N" 'weights 8 1.0 0.7 1.0 0.8 1.3 0.9 1.0 0.6'      "$COUNT"
expect_error nurbs-zero       iga-nurbs-weights.in "$N" 'weights 9 1.0 0.7 1.0 0.8 0.0 0.9 1.0 0.6 1.2'  "$POSITIVE"
expect_error nurbs-negative   iga-nurbs-weights.in "$N" 'weights 9 1.0 0.7 1.0 0.8 -1.3 0.9 1.0 0.6 1.2' "$POSITIVE"
expect_error nurbs-nan        iga-nurbs-weights.in "$N" 'weights 9 1.0 0.7 1.0 0.8 nan 0.9 1.0 0.6 1.2'  "$POSITIVE"
expect_error tspline-count    iga-tspline-nip.in   "$T" 'weights 15 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1'       "$COUNT"
expect_error tspline-negative iga-tspline-nip.in   "$T" 'weights 16 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 -1'    "$POSITIVE"

[ "$failures" -eq 0 ]
