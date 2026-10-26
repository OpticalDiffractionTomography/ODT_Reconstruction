#!/bin/bash
# test_tomo_process.sh — tests for the tomo_process CLI: --restart, --resume,
# --cancel and argument dispatch. The orchestrator is replaced by a fake that
# just records that it started; scancel/sbatch are mocked.
#
#   bash tests/test_tomo_process.sh

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
SRC_FILE="${REPO_ROOT}/tomo_process"

BG_PIDS=()
trap 'kill "${BG_PIDS[@]}" 2>/dev/null; rm -rf "${TMP_ROOT}"' EXIT

# A run dir with saved parameters, plus a fake repo whose orchestrator.sh only
# records its start and then idles (so --restart has something to start).
prepare() {
    # die/info are one-liners in tomo_process (not extractable); same behaviour
    die()  { echo "ERROR: $*" >&2; exit 1; }
    info() { echo "[tomo_process] $*"; }
    load_fns cmd_resume cmd_restart cmd_cancel
    SCRATCH_ROOT="${TEST_DIR}/scratch"
    RUN="${SCRATCH_ROOT}/run1"
    mkdir -p "${RUN}/state" "${RUN}/logs"
    PIDFILE="${RUN}/state/orchestrator.pid"
    printf 'c1\t/d\tf.mat\n' > "${RUN}/state/chunks.txt"
    SCRIPT_DIR="${TEST_DIR}/fakerepo"
    mkdir -p "${SCRIPT_DIR}/scripts"
    cat > "${SCRIPT_DIR}/scripts/orchestrator.sh" <<'EOF'
#!/bin/bash
echo "started pid=$$ run=${TOMO_RUN_ID} version=${TOMO_VERSION}" >> "${TOMO_SCRATCH_ROOT}/${TOMO_RUN_ID}/started.log"
sleep 20
EOF
    cat > "${RUN}/state/run_env.sh" <<EOF
TOMO_RUN_ID=run1; TOMO_INPUT_PATH=p; TOMO_EMAIL=me@example.org; TOMO_MAX_JOBS=2
TOMO_MINS_PER_SAMPLE=5; TOMO_REPO_DIR=${SCRIPT_DIR}; TOMO_SCRATCH_ROOT=${SCRATCH_ROOT}
TOMO_DATA_MOUNT=/d; TOMO_RESULTS_MOUNT=/r; TOMO_POLL_INTERVAL=1
TOMO_SLURM_PARTITION=p; TOMO_SLURM_GRES=g; TOMO_SLURM_CPUS=1; TOMO_SLURM_MEM=1G
TOMO_SLURM_TIME=01:00:00; TOMO_NM=1.337; TOMO_VERSION=alice
EOF
}
wait_for_file() {  # path
    local i
    for i in $(seq 1 50); do [[ -s "$1" ]] && return 0; sleep 0.1; done
    return 1
}
alive() { kill -0 "$1" 2>/dev/null; }
start_dummy_orchestrator() {  # sets OLD_PID and writes it to the pid file
    sleep 20 &
    OLD_PID=$!; BG_PIDS+=("${OLD_PID}")
    echo "${OLD_PID}" > "${PIDFILE}"
}
new_pid_started() { wait_for_file "${RUN}/started.log" && BG_PIDS+=("$(cat "${PIDFILE}")"); }

test_restart_stops_old_orchestrator_and_starts_a_new_one() {
    prepare; start_dummy_orchestrator
    cmd_restart run1 > "${TEST_DIR}/out" 2>&1
    ! alive "${OLD_PID}"
    new_pid_started
    local new; new="$(cat "${PIDFILE}")"
    [[ "${new}" != "${OLD_PID}" ]]
    alive "${new}"
    assert_contains "${RUN}/started.log" "run=run1 version=alice"   # saved params are reused
    kill "${new}" 2>/dev/null || true
}
test_restart_works_when_pid_file_is_missing() {
    prepare
    cmd_restart run1 > "${TEST_DIR}/out" 2>&1
    new_pid_started
    kill "$(cat "${PIDFILE}")" 2>/dev/null || true
}
test_restart_works_when_recorded_pid_is_stale() {
    prepare
    echo 999999 > "${PIDFILE}"
    cmd_restart run1 > "${TEST_DIR}/out" 2>&1
    new_pid_started
    kill "$(cat "${PIDFILE}")" 2>/dev/null || true
}
test_restart_of_unknown_run_fails() {
    prepare
    ( cmd_restart nope ) > "${TEST_DIR}/out" 2>&1 && return 1
    assert_contains "${TEST_DIR}/out" "Run not found: nope"
}
test_resume_refuses_while_orchestrator_is_running() {
    prepare; start_dummy_orchestrator
    ( cmd_resume run1 ) > "${TEST_DIR}/out" 2>&1 && return 1
    assert_contains "${TEST_DIR}/out" "already running"
    alive "${OLD_PID}"        # and it must not have been touched
    kill "${OLD_PID}" 2>/dev/null || true
}
test_resume_requires_saved_parameters() {
    prepare
    rm "${RUN}/state/run_env.sh"
    ( cmd_resume run1 ) > "${TEST_DIR}/out" 2>&1 && return 1
    assert_contains "${TEST_DIR}/out" "run_env.sh"
}
test_cancel_stops_orchestrator_and_cancels_notify_job_before_chunk_jobs() {
    prepare; start_dummy_orchestrator
    echo 7000 > "${RUN}/state/start_notify_jid"
    printf 'c1\t7001\nc2\t7002\n' > "${RUN}/state/active_jobs.txt"
    SCANCELLED="${TEST_DIR}/scancelled"; : > "${SCANCELLED}"
    scancel() { echo "$1" >> "${SCANCELLED}"; }
    sbatch() { return 0; }
    cmd_cancel run1 > "${TEST_DIR}/out" 2>&1
    ! alive "${OLD_PID}"
    [[ -f "${RUN}/state/cancelled" ]]
    # the "Started" notification must go first, or cancelling chunks would release it
    assert_eq "7000 7001 7002" "$(tr '\n' ' ' < "${SCANCELLED}" | sed 's/ $//')"
}

# ── CLI dispatch (runs the real script) ───────────────────────────────────────
test_help_lists_restart() {
    bash "${REPO_ROOT}/tomo_process" --help > "${TEST_DIR}/out" 2>&1
    assert_contains "${TEST_DIR}/out" "--restart <run_id>"
}
test_restart_without_run_id_is_an_error() {
    bash "${REPO_ROOT}/tomo_process" --restart > "${TEST_DIR}/out" 2>&1 && return 1
    assert_contains "${TEST_DIR}/out" "--restart requires a run_id"
}
test_no_arguments_prints_usage_and_fails() {
    bash "${REPO_ROOT}/tomo_process" > "${TEST_DIR}/out" 2>&1 && return 1
    assert_contains "${TEST_DIR}/out" "Usage:"
}

run_all "tomo_process"
