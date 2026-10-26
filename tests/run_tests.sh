#!/bin/bash
# run_tests.sh — run all tests. Run from anywhere.
#   bash tests/run_tests.sh
#
# Shell tests always run. MATLAB unit tests run when `matlab` is on PATH
# (e.g. after `module load matlab/R2023a` on the cluster); otherwise skipped.

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
status=0

echo "== shell tests =="
bash "${ROOT}/tests/test_orchestrator.sh" || status=1

echo
echo "== MATLAB tests =="
if command -v matlab > /dev/null 2>&1; then
    matlab -batch "cd('${ROOT}'); r = runtests('tests'); disp(table(r)); exit(any([r.Failed] | [r.Incomplete]))" \
        || status=1
else
    echo "matlab not found on PATH — skipped"
fi

exit "${status}"
