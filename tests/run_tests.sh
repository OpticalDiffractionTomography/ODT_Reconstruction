#!/bin/bash
# run_tests.sh — run all tests. Run from anywhere.
#   bash tests/run_tests.sh
#
# Every tests/test_*.sh runs. MATLAB unit tests (tests/t*.m) run when `matlab`
# is on PATH (e.g. after `module load matlab/R2023a` on the cluster); otherwise
# they are skipped.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
status=0

echo "== shell tests =="
for t in "${ROOT}"/tests/test_*.sh; do
    echo "-- $(basename "${t}")"
    bash "${t}" || status=1
    echo
done

echo "== MATLAB tests =="
if command -v matlab > /dev/null 2>&1; then
    matlab -batch "cd('${ROOT}'); r = runtests('tests'); disp(table(r)); exit(any([r.Failed] | [r.Incomplete]))" \
        || status=1
else
    echo "matlab not found on PATH — skipped"
fi

exit "${status}"
