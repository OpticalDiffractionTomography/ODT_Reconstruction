#!/bin/bash
# test_job_template.sh — runs the real scripts/bash_template.sh end to end with a
# fake `matlab` (and `module`) on PATH, to check the DONE/FAIL markers and
# stage sequencing without MATLAB or SLURM.
#
#   bash tests/test_job_template.sh
#
# The fake matlab is steered per stage by FAKE_S1 / FAKE_S2:
#   ok            exit 0 (stage 1 also creates a Field_a.mat)
#   fail          exit 1, no outputs
#   fail_partial  stage 1 only: create Field_a.mat, then exit 1 (a crash after
#                 one sample was finished)
#   slow_ok       like ok, but takes ~2 s and ignores SIGTERM (a long stage)

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"

# Build a job script from the template, a fake bin dir, and the env to run it
prepare() {
    DATA="${TEST_DIR}/chunk/data"
    mkdir -p "${DATA}/field_retrieval" "${TEST_DIR}/bin"
    DONE="${TEST_DIR}/chunk/DONE"; FAILM="${TEST_DIR}/chunk/FAIL"
    FAKE_LOG="${TEST_DIR}/matlab_calls.log"; FAKE_CMDS="${TEST_DIR}/matlab_cmds.log"
    OUT="${TEST_DIR}/job.out"
    : > "${FAKE_LOG}"; : > "${FAKE_CMDS}"
    sed -e "s|__JOB_NAME__|job|g" -e "s|__LOG_DIR__|${TEST_DIR}|g" \
        -e "s|__DATA_DIR__|${DATA}|g" -e "s|__CHUNK_DONE__|${DONE}|g" \
        -e "s|__CHUNK_FAIL__|${FAILM}|g" -e "s|__REPO_DIR__|${REPO_ROOT}|g" \
        -e "s|__PARTITION__|p|g" -e "s|__GRES__|g|g" -e "s|__CPUS__|1|g" \
        -e "s|__MEM__|1G|g" -e "s|__TIME__|01:00:00|g" -e "s|__NM__|1.337|g" \
        -e "s|__VERSION__|${VERSION:-default}|g" -e "s|__BATCH_NAME__|batchX|g" \
        -e "s|__TIFF_OFFSET__|${OFFSET:-0}|g" \
        "${REPO_ROOT}/scripts/bash_template.sh" > "${TEST_DIR}/job.sh"
    printf '#!/bin/bash\nexit 0\n' > "${TEST_DIR}/bin/module"
    cat > "${TEST_DIR}/bin/matlab" <<'EOF'
#!/bin/bash
cmd="$2"
case "$cmd" in
    *field_Retrieval.m*)         stage=1 ;;
    *tomogram_Reconstruction.m*) stage=2 ;;
    *) echo "fake matlab: unknown stage" >&2; exit 98 ;;
esac
echo "${stage}" >> "${FAKE_LOG}"
echo "${cmd}" >> "${FAKE_CMDS}"
mode_var="FAKE_S${stage}"; mode="${!mode_var:-ok}"
make_field() { touch "${FAKE_DATA}/field_retrieval/Field_a.mat"; }
case "${mode}" in
    ok)           [[ "${stage}" -eq 1 ]] && make_field; exit 0 ;;
    fail)         exit 1 ;;
    fail_partial) make_field; exit 1 ;;
    slow_ok)      trap '' TERM; sleep 2; [[ "${stage}" -eq 1 ]] && make_field; exit 0 ;;
esac
EOF
    chmod +x "${TEST_DIR}/bin/matlab" "${TEST_DIR}/bin/module"
    export PATH="${TEST_DIR}/bin:${PATH}" FAKE_LOG FAKE_CMDS FAKE_DATA="${DATA}"
}
run_job() { bash "${TEST_DIR}/job.sh" > "${OUT}" 2>&1; }
stages_run() { tr '\n' ' ' < "${FAKE_LOG}" | sed 's/ $//'; }

test_success_writes_done_marker_and_runs_both_stages() {
    prepare
    run_job
    [[ -f "${DONE}" && ! -f "${FAILM}" ]]
    assert_eq "1 2" "$(stages_run)"
}
test_stage1_crash_after_partial_output_still_runs_stage2_then_fails() {
    prepare; export FAKE_S1=fail_partial
    ! run_job
    assert_eq "1 2" "$(stages_run)" "stage 2 must still reconstruct the finished Field files"
    [[ -f "${FAILM}" && ! -f "${DONE}" ]]
}
test_stage1_failure_without_output_skips_stage2_and_fails() {
    prepare; export FAKE_S1=fail
    ! run_job
    assert_eq "1" "$(stages_run)"
    assert_contains "${OUT}" "skipping tomogram_Reconstruction"
    [[ -f "${FAILM}" && ! -f "${DONE}" ]]
}
test_stage2_failure_fails_the_chunk() {
    prepare; export FAKE_S2=fail
    ! run_job
    assert_eq "1 2" "$(stages_run)"
    [[ -f "${FAILM}" && ! -f "${DONE}" ]]
}
test_both_stages_failing_is_still_reported_once_as_failed() {
    prepare; export FAKE_S1=fail_partial FAKE_S2=fail
    ! run_job
    [[ -f "${FAILM}" && ! -f "${DONE}" ]]
}
test_job_exit_code_is_nonzero_on_failure_and_zero_on_success() {
    prepare
    run_job
    prepare; export FAKE_S2=fail
    ! run_job
}
test_stage_arguments_reach_matlab() {
    VERSION=alice OFFSET=3 prepare
    run_job
    assert_contains "${FAKE_CMDS}" "spath='${DATA}'"
    assert_contains "${FAKE_CMDS}" "pipeline_version='alice'"
    assert_contains "${FAKE_CMDS}" "n_m=1.337"
    assert_contains "${FAKE_CMDS}" "batch_name='batchX'"
    assert_contains "${FAKE_CMDS}" "tiff_offset=3"
}
test_sigterm_writes_fail_marker_even_if_stages_would_succeed() {
    # SLURM sends SIGTERM at the walltime limit; the job must leave a FAIL
    # marker instead of vanishing with neither marker.
    prepare; export FAKE_S1=slow_ok FAKE_S2=slow_ok
    bash "${TEST_DIR}/job.sh" > "${OUT}" 2>&1 &
    local pid=$! i
    for i in $(seq 1 50); do [[ -s "${FAKE_LOG}" ]] && break; sleep 0.1; done
    [[ -s "${FAKE_LOG}" ]]          # matlab is running
    kill -TERM "${pid}"
    wait "${pid}" || true
    [[ -f "${FAILM}" && ! -f "${DONE}" ]]
}

run_all "bash_template.sh"
