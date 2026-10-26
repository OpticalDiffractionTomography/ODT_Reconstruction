#!/bin/bash
# lib.sh — tiny test harness shared by tests/test_*.sh. Source it, set SRC_FILE
# to the script under test, define test_* functions, then call run_all.
#
# Scripts such as orchestrator.sh are long-running and cannot be sourced, so
# tests pull the functions they need out by name (load_fns) and run each test
# in a subshell with a fresh scratch directory (setup_env).

set -u
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_ROOT="$(mktemp -d)"
trap 'rm -rf "${TMP_ROOT}"' EXIT
PASS=0
FAIL=0
SRC_FILE=""

# Print the source of function $1 from $SRC_FILE
extract_fn() {
    awk -v n="$1" 'index($0, n "() {") == 1 {p=1} p {print} p && $0 == "}" {exit}' "${SRC_FILE}"
}

# Load the named functions into the current shell
load_fns() {
    local f src
    for f in "$@"; do
        src="$(extract_fn "${f}")"
        [[ -n "${src}" ]] || { echo "cannot extract function ${f} from ${SRC_FILE}"; exit 99; }
        eval "${src}"
    done
}

# Fresh orchestrator-style state for one test: RUN_SCRATCH, STATE_DIR, the
# state files, and mount dirs under a new temp dir.
setup_env() {
    local d
    d="$(mktemp -d -p "${TMP_ROOT}")"
    TEST_DIR="${d}"
    RUN_SCRATCH="${d}/run"
    STATE_DIR="${RUN_SCRATCH}/state"
    LOG_DIR="${RUN_SCRATCH}/logs"
    JOBS_DIR="${RUN_SCRATCH}/jobs"
    mkdir -p "${STATE_DIR}" "${LOG_DIR}" "${JOBS_DIR}"
    CHUNKS_FILE="${STATE_DIR}/chunks.txt"
    ACTIVE_FILE="${STATE_DIR}/active_jobs.txt"
    DONE_FILE="${STATE_DIR}/done_chunks.txt"
    FAIL_FILE="${STATE_DIR}/failed_chunks.txt"
    COPYFAIL_FILE="${STATE_DIR}/copyfail_chunks.txt"
    SKIPPED_FILE="${STATE_DIR}/skipped_files.txt"
    RETRIED_FILE="${STATE_DIR}/retried_chunks.txt"
    SUBMIT_TIMES_FILE="${STATE_DIR}/submit_times.txt"
    DURATIONS_FILE="${STATE_DIR}/chunk_durations.txt"
    TOMO_DATA_MOUNT="${d}/data"
    TOMO_RESULTS_MOUNT="${d}/results"
    TOMO_RUN_ID="run1"
    : > "${CHUNKS_FILE}"; : > "${ACTIVE_FILE}"; : > "${DONE_FILE}"; : > "${FAIL_FILE}"
    : > "${COPYFAIL_FILE}"; : > "${SKIPPED_FILE}"; : > "${RETRIED_FILE}"
    : > "${SUBMIT_TIMES_FILE}"; : > "${DURATIONS_FILE}"
    log() { :; }
}

# run NAME FUNCTION — run a test in a subshell and tally the result
run() {
    local name="$1" fn="$2" out rc
    # Capture the status rather than testing the subshell in an `if`: bash
    # ignores `set -e` inside a command that is part of a condition, which
    # would let every assertion but the last pass silently.
    out="$( ( set -e; setup_env; "${fn}" ) 2>&1 )"; rc=$?
    if [[ "${rc}" -eq 0 ]]; then
        PASS=$(( PASS + 1 )); echo "  ok    ${name}"
    else
        FAIL=$(( FAIL + 1 )); echo "  FAIL  ${name}"; echo "${out}" | sed 's/^/        /'
    fi
}

assert_eq() {  # expected actual [message]
    [[ "$1" == "$2" ]] || { echo "expected: [$1]"; echo "actual:   [$2]"; echo "${3:-assertion failed}"; return 1; }
}
assert_file_lines() {  # file expected_count
    local n; n=$(wc -l < "$1")
    [[ "${n}" -eq "$2" ]] || { echo "$1 has ${n} line(s), expected $2"; cat "$1"; return 1; }
}
assert_contains() {  # file pattern [message]
    grep -qE -- "$2" "$1" || { echo "${3:-$1 does not contain /$2/}"; echo "--- $1:"; cat "$1"; return 1; }
}
assert_not_contains() {  # file pattern [message]
    ! grep -qE -- "$2" "$1" || { echo "${3:-$1 unexpectedly contains /$2/}"; echo "--- $1:"; cat "$1"; return 1; }
}

# run_all TITLE — run every test_* function defined in the calling script
run_all() {
    echo "$1"
    local fn
    for fn in $(declare -F | awk '{print $3}' | grep '^test_'); do
        run "${fn#test_}" "${fn}"
    done
    echo
    echo "${PASS} passed, ${FAIL} failed"
    [[ "${FAIL}" -eq 0 ]]
}
