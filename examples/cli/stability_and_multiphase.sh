#!/usr/bin/env bash
# The ADR-0033 CLI surface: stability, PC-SAFT, and more than two phases.
# Run from the repo root:  bash examples/cli/stability_and_multiphase.sh
# PYTHON selects the interpreter (default: python3); `chemthermo ...` is the
# same as `$PYTHON -m chemthermo ...` once the package is installed.
set -euo pipefail
PYTHON="${PYTHON:-python3}"
ct() { "$PYTHON" -m chemthermo "$@"; }

echo "== 1. A stable feed (Peng-Robinson). 'stable' is bounded by the trial set."
ct stability-tp --components Methane,Ethane --z 0.5,0.5 --temperature-k 300 --pressure-pa 100000

echo
echo "== 2. An unstable feed: the incipient vapour it would form."
ct stability-tp --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000

echo
echo "== 3. The same feed with PC-SAFT, as JSON (verdict and scope only)."
ct stability-tp --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000 \
    --eos pc-saft --format json | "$PYTHON" -c \
    'import json, sys; r = json.load(sys.stdin)["result"]; print(r["status"], r["stability_scope"], r["tpd_min"])'

echo
echo "== 4. Its PC-SAFT flash."
ct tp-flash --components Methane,n-Hexane --z 0.5,0.5 --temperature-k 300 --pressure-pa 2000000 --eos pc-saft

echo
echo "== 5. Three liquids: Peng-Robinson water / ethanol / n-hexane at 280 K, 1 atm."
ct tp-flash --components Water,Ethanol,n-Hexane --z 0.2,0.4,0.4 --temperature-k 280 --pressure-pa 101325 \
    --max-phases 3

echo
echo "== 6. The same state with a two-phase budget is refused (exit code 3), not faked."
set +e
ct tp-flash --components Water,Ethanol,n-Hexane --z 0.2,0.4,0.4 --temperature-k 280 --pressure-pa 101325 \
    --max-phases 2 2>/dev/null
code=$?
set -e
echo "exit code: $code"
test "$code" -eq 3
