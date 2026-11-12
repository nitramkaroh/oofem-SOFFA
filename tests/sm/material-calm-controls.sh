# Material CALM: the full norm and HPC agree when they select the same dofs,
# and a failed Newton solve must be reported as a failed engineering-model step.
# Run by ctest as: bash material-calm-controls.sh <oofem executable>
OOFEM=$1
case "$OOFEM" in
    /* | ?:*) ;;
    *) OOFEM="$PWD/$OOFEM" ;;
esac
INPUT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) || exit 1
INPUT="$INPUT_DIR/material-calm-psi.in"
MIXED_INPUT="$INPUT_DIR/surface-arclength-mixed.in"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/oofem-material-calm.XXXXXXXX") || exit 1
# These cases write only flat input, output and log files; no recursive cleanup is needed.
cleanup() {
    cd / || return
    rm -f "$WORK"/*.in "$WORK"/*.out "$WORK"/*.log
    rmdir "$WORK"
}
trap cleanup EXIT
cd "$WORK" || exit 1

failures=0

# The same analytic checks as the HPC input must pass with the full displacement norm.
sed -e '1s/.*/full-psi.out/' -e 's/ hpc 2 2 1//' "$INPUT" > full-psi.in
if "$OOFEM" -f full-psi.in > full-psi.log 2>&1; then
    if grep -q 'Checking rules' full-psi.log && ! grep -qi 'nan' full-psi.out; then
        echo 'OK   full-norm-psi: analytic equilibrium and arc length'
    else
        echo 'FAIL full-norm-psi: checks were not evaluated or the result contains NaN'
        cat full-psi.log
        failures=$((failures + 1))
    fi
else
    echo 'FAIL full-norm-psi: the analytic control case did not finish'
    cat full-psi.log
    failures=$((failures + 1))
fi

# Convergence attained on the last allowed iteration is success, not an exhausted iteration budget.
# Psi=0 fixes u=-0.02 in this one-dof example, and the first correction finds Lambda exactly.
sed -e '1s/.*/last-iteration.out/' -e 's/MaxIter 30/MaxIter 1 maxrestarts 0/' \
    "$MIXED_INPUT" > last-iteration.in
if "$OOFEM" -f last-iteration.in > last-iteration.log 2>&1; then
    if grep -q 'Checking rules' last-iteration.log &&
       grep -q 'Total 0 error(s) and 0 warning(s)' last-iteration.log; then
        echo 'OK   last-iteration: convergence on MaxIter is accepted without a warning'
    else
        echo 'FAIL last-iteration: a converged final iteration was treated as failure'
        cat last-iteration.log
        failures=$((failures + 1))
    fi
else
    echo 'FAIL last-iteration: a converged final iteration stopped the analysis'
    cat last-iteration.log
    failures=$((failures + 1))
fi

# A step of 0.2 needs more than three iterations. After one restart with the
# reduction factor 0.25, the accepted radius is 0.05. Independently solving
# u^2 + (s - 1/s)^2 = 0.05^2 and Lambda = 1 - 2(s - 1/s), s = 1 + u,
# gives the equilibrium below. This checks that a retry restores its starting
# load level, displacements and material state rather than keeping the failed trial.
awk '
    /^#%BEGIN_CHECK%/ { exit }
    NR == 1 { print "restarted.out"; next }
    {
        gsub(" hpc 2 2 1", "")
        gsub("stepLength 0.02", "stepLength 0.2")
        gsub("MaxIter 30", "MaxIter 3 maxrestarts 4")
        print
    }
' "$INPUT" > restarted.in
cat >> restarted.in <<'EOF'
#%BEGIN_CHECK% tolerance 1.e-10
#NODE tStep 1 number 2 dof 1 unknown d value -0.022159584152275758
#LOADLEVEL tStep 1 value 1.0896426869342548
#%END_CHECK%
EOF
if "$OOFEM" -f restarted.in > restarted.log 2>&1; then
    if grep -q 'Iteration Reset' restarted.log &&
       grep -q 'Checking rules' restarted.log &&
       grep -q 'Total 0 error(s) and 0 warning(s)' restarted.log; then
        echo 'OK   restarted: reduced step reaches its analytic equilibrium'
    else
        echo 'FAIL restarted: the step did not restart and converge cleanly'
        cat restarted.log
        failures=$((failures + 1))
    fi
else
    echo 'FAIL restarted: the reduced step did not converge'
    cat restarted.log
    failures=$((failures + 1))
fi

# A larger full-norm step with only one allowed Newton iteration cannot meet 1e-11.
# No errorcheck module is used: this specifically checks propagation of CR_DIVERGED_ITS
# from the material CALM solver to StaticStructural, rather than an unrelated golden-value error.
awk '
    /^#%BEGIN_CHECK%/ { exit }
    /^errorcheck[[:space:]]*$/ { next }
    NR == 1 { print "nonconverged.out"; next }
    {
        gsub(" hpc 2 2 1", "")
        gsub("stepLength 0.02", "stepLength 0.2")
        gsub("MaxIter 30", "MaxIter 1 maxrestarts 0")
        gsub("nmodules 1", "nmodules 0")
        print
    }
' "$INPUT" > nonconverged.in
if "$OOFEM" -f nonconverged.in > nonconverged.log 2>&1; then
    echo 'FAIL nonconverged: the analysis accepted an unconverged step'
    cat nonconverged.log
    failures=$((failures + 1))
elif grep -q 'Nonlinear solver failed to converge at step' nonconverged.log; then
    echo 'OK   nonconverged: failed solver status reached the engineering model'
else
    echo 'FAIL nonconverged: the analysis stopped for another reason'
    cat nonconverged.log
    failures=$((failures + 1))
fi

# A selected parameter with zero amplitude provides no continuation direction, even though
# the other membrane still supplies a finite stiffness. Stop before dividing by its zero norm.
awk '
    /^#%BEGIN_CHECK%/ { exit }
    /^errorcheck[[:space:]]*$/ { next }
    NR == 1 { print "zero-direction.out"; next }
    {
        gsub("gamma1 0.5", "gamma1 0.0")
        gsub("nmodules 1", "nmodules 0")
        print
    }
' "$INPUT" > zero-direction.in
if "$OOFEM" -f zero-direction.in > zero-direction.log 2>&1; then
    echo 'FAIL zero-direction: the analysis accepted a zero continuation direction'
    cat zero-direction.log
    failures=$((failures + 1))
elif grep -q 'no nonzero direction' zero-direction.log; then
    echo 'OK   zero-direction: the zero continuation direction was rejected'
else
    echo 'FAIL zero-direction: the analysis stopped for another reason'
    cat zero-direction.log
    failures=$((failures + 1))
fi

[ "$failures" -eq 0 ]
