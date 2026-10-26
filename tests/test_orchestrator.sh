#!/bin/bash
# test_orchestrator.sh — unit tests for the failure-handling logic in
# scripts/orchestrator.sh. Needs only bash + awk; squeue/rsync are mocked.
#
#   bash tests/test_orchestrator.sh
#
# The orchestrator is a long-running script and cannot be sourced, so each
# test extracts the functions it needs by name and runs them in a subshell
# with a fresh scratch directory.

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
SRC_FILE="${REPO_ROOT}/scripts/orchestrator.sh"

# ── job_active ────────────────────────────────────────────────────────────────
test_job_active_running() {
    load_fns job_active
    squeue() { echo "12345 RUNNING"; }
    job_active 12345
}
test_job_active_finished_empty_output() {
    load_fns job_active
    squeue() { return 0; }   # squeue exits 0 with no output for a just-finished job
    ! job_active 12345
}
test_job_active_purged_invalid_id() {
    load_fns job_active
    squeue() { echo "slurm_load_jobs error: Invalid job id specified"; return 1; }
    ! job_active 12345
}
test_job_active_slurm_outage_counts_as_active() {
    load_fns job_active
    squeue() { echo "slurm_load_jobs error: Unable to contact slurm controller"; return 1; }
    job_active 12345
}

# ── sample_index ──────────────────────────────────────────────────────────────
test_sample_index_positions_per_dir() {
    load_fns sample_index
    printf '%s\n' /m/A/s001.mat /m/A/s002.mat /m/A/sub/x.mat /m/A/s003.mat /m/B/s001.mat \
        > "${STATE_DIR}/pending_files.txt"
    assert_eq 1 "$(sample_index /m/A /m/A/s001.mat)"
    assert_eq 3 "$(sample_index /m/A /m/A/s003.mat)" "files in subdirs must not be counted"
    assert_eq 1 "$(sample_index /m/B /m/B/s001.mat)" "numbering restarts per source dir"
    assert_eq "" "$(sample_index /m/A /m/A/nope.mat)"
}

# ── unfinished_files ──────────────────────────────────────────────────────────
test_unfinished_files_recognises_both_naming_schemes() {
    load_fns unfinished_files
    local rs="${RUN_SCRATCH}/c1/data/field_retrieval"
    mkdir -p "${rs}"
    touch "${rs}/Tomogram_Field_s1_Tomog_bg_1.mat" "${rs}/Tomogram_Field_s2_Tomog.mat"
    local out
    out="$(unfinished_files c1 "/x/s1_Tomog.mat,/x/s2_Tomog.mat,/x/s3_Tomog.mat")"
    assert_eq "/x/s3_Tomog.mat" "${out}"
}
test_unfinished_files_ignores_part_files() {
    load_fns unfinished_files
    local rs="${RUN_SCRATCH}/c1/data/field_retrieval"
    mkdir -p "${rs}"
    touch "${rs}/Tomogram_Field_s1_Tomog.mat.part"   # incomplete atomic save
    assert_eq "/x/s1_Tomog.mat" "$(unfinished_files c1 "/x/s1_Tomog.mat")"
}
test_unfinished_files_excludes_reported_skips() {
    load_fns unfinished_files
    mkdir -p "${RUN_SCRATCH}/c1/data/field_retrieval"
    printf 'c1\ts2_Tomog.mat\tcorrupt\n' > "${SKIPPED_FILE}"
    assert_eq "/x/s1_Tomog.mat,/x/s3_Tomog.mat" \
        "$(unfinished_files c1 "/x/s1_Tomog.mat,/x/s2_Tomog.mat,/x/s3_Tomog.mat")"
}
test_unfinished_files_skip_of_other_chunk_does_not_hide_file() {
    load_fns unfinished_files
    mkdir -p "${RUN_SCRATCH}/c1/data/field_retrieval"
    printf 'c9\ts1_Tomog.mat\tcorrupt\n' > "${SKIPPED_FILE}"
    assert_eq "/x/s1_Tomog.mat" "$(unfinished_files c1 "/x/s1_Tomog.mat")"
}

# ── requeue_unfinished ────────────────────────────────────────────────────────
test_requeue_nothing_remaining() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A ""
    assert_eq none "${REQUEUE_RESULT}"
    assert_file_lines "${CHUNKS_FILE}" 0
}
test_requeue_bisects_even_count() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A "/m/A/a,/m/A/b,/m/A/c,/m/A/d"
    assert_eq retry "${REQUEUE_RESULT}"
    assert_file_lines "${CHUNKS_FILE}" 2
    assert_eq "/m/A/a,/m/A/b" "$(sed -n 1p "${CHUNKS_FILE}" | cut -f3)"
    assert_eq "/m/A/c,/m/A/d" "$(sed -n 2p "${CHUNKS_FILE}" | cut -f3)"
}
test_requeue_bisects_odd_count_without_losing_files() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A "/m/A/a,/m/A/b,/m/A/c"
    assert_eq "/m/A/a,/m/A/b" "$(sed -n 1p "${CHUNKS_FILE}" | cut -f3)"
    assert_eq "/m/A/c" "$(sed -n 2p "${CHUNKS_FILE}" | cut -f3)"
}
test_requeue_new_chunks_get_unique_ids_and_source_dir() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A "/m/A/a,/m/A/b"
    requeue_unfinished c2 /m/B "/m/B/a,/m/B/b"
    local ids; ids="$(cut -f1 "${CHUNKS_FILE}" | sort -u | wc -l)"
    assert_eq 4 "${ids}" "retry chunk ids must be unique"
    assert_eq /m/B "$(sed -n 3p "${CHUNKS_FILE}" | cut -f2)"
}
test_requeue_single_file_gets_one_extra_attempt_then_is_skipped() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A "/m/A/bad.mat"
    assert_eq retry "${REQUEUE_RESULT}" "first failure should be retried"
    assert_file_lines "${CHUNKS_FILE}" 1
    assert_file_lines "${SKIPPED_FILE}" 0
    requeue_unfinished retry0001 /m/A "/m/A/bad.mat"
    assert_eq skipped "${REQUEUE_RESULT}" "second failure should give up"
    assert_file_lines "${CHUNKS_FILE}" 1 "no further chunk may be queued"
    assert_file_lines "${SKIPPED_FILE}" 1
    grep -q "bad.mat" "${SKIPPED_FILE}"
}
test_requeue_records_retries_for_summary() {
    load_fns requeue_unfinished
    requeue_unfinished c1 /m/A "/m/A/a,/m/A/b"
    assert_file_lines "${RETRIED_FILE}" 1
}

# ── copy_chunk_results ────────────────────────────────────────────────────────
make_chunk_results() {  # cid
    local cid="$1"
    mkdir -p "${RUN_SCRATCH}/${cid}/data/field_retrieval"
    printf '%s\t%s\t%s\n' /dummy /dummy /dummy > "${CHUNKS_FILE}"
    printf '%s\t%s\t%s\n' "${cid}" "${TOMO_DATA_MOUNT}/proj/exp" "f.mat" > "${CHUNKS_FILE}"
}
test_copy_empty_results_is_not_an_error() {
    load_fns copy_chunk_results
    make_chunk_results c1
    rsync() { echo "rsync must not be called" >&2; return 1; }
    copy_chunk_results c1
}
test_copy_suffixes_logs_and_reports_skipped_files() {
    load_fns copy_chunk_results
    make_chunk_results c1
    local rs="${RUN_SCRATCH}/c1/data/field_retrieval"
    touch "${rs}/Tomogram_Field_s1.mat" "${rs}/field_retrieval.log"
    printf 's2_Tomog.mat\tfield_Retrieval\tcannot read MAT-file\n' > "${rs}/skipped_files.txt"
    rsync() { return 0; }
    copy_chunk_results c1
    [[ -f "${rs}/field_retrieval_c1.log" && -f "${rs}/skipped_files_c1.txt" ]]
    assert_file_lines "${SKIPPED_FILE}" 1
    assert_eq "c1" "$(cut -f1 "${SKIPPED_FILE}")"
    assert_eq "s2_Tomog.mat" "$(cut -f2 "${SKIPPED_FILE}")"
}
test_copy_excludes_part_files_and_overwrites() {
    load_fns copy_chunk_results
    make_chunk_results c1
    touch "${RUN_SCRATCH}/c1/data/field_retrieval/a.mat"
    rsync() { echo "$*" > "${STATE_DIR}/rsync_args"; return 0; }
    copy_chunk_results c1
    grep -q -- "--exclude=\*.part" "${STATE_DIR}/rsync_args"
    grep -q -- "--ignore-times" "${STATE_DIR}/rsync_args"
}
test_copy_failure_returns_nonzero_and_keeps_results() {
    load_fns copy_chunk_results
    make_chunk_results c1
    touch "${RUN_SCRATCH}/c1/data/field_retrieval/a.mat"
    rsync() { return 23; }
    ! copy_chunk_results c1
    [[ -f "${RUN_SCRATCH}/c1/data/field_retrieval/a.mat" ]]
}

# ── collect_finished ──────────────────────────────────────────────────────────
load_collect_fns() {
    load_fns fmt_duration record_chunk_duration job_active unfinished_files \
        requeue_unfinished copy_chunk_results collect_finished
    rsync() { return 0; }
    squeue() { return 0; }
}
make_active_chunk() {  # cid jid files_csv
    mkdir -p "${RUN_SCRATCH}/$1/data/field_retrieval"
    printf '%s\t%s\t%s\n' "$1" "${TOMO_DATA_MOUNT}/proj/exp" "$3" >> "${CHUNKS_FILE}"
    printf '%s\t%s\n' "$1" "$2" >> "${ACTIVE_FILE}"
}
test_collect_done_chunk_is_copied_and_removed() {
    load_collect_fns
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    touch "${RUN_SCRATCH}/c1/DONE" "${RUN_SCRATCH}/c1/data/field_retrieval/Tomogram_Field_s1_Tomog.mat"
    collect_finished
    assert_file_lines "${DONE_FILE}" 1
    [[ ! -d "${RUN_SCRATCH}/c1" ]]
}
test_collect_failed_chunk_salvages_and_requeues_unfinished() {
    load_collect_fns
    make_active_chunk c1 100 "/x/s1_Tomog.mat,/x/s2_Tomog.mat,/x/s3_Tomog.mat,/x/s4_Tomog.mat"
    local rs="${RUN_SCRATCH}/c1/data/field_retrieval"
    touch "${RUN_SCRATCH}/c1/FAIL" "${rs}/Tomogram_Field_s1_Tomog.mat" "${rs}/Tomogram_Field_s2_Tomog.mat"
    collect_finished
    # s1/s2 finished -> not re-queued; s3/s4 re-queued; parent not counted as lost
    assert_file_lines "${CHUNKS_FILE}" 3   # original + 2 retry halves
    assert_eq "/x/s3_Tomog.mat" "$(sed -n 2p "${CHUNKS_FILE}" | cut -f3)"
    assert_eq "/x/s4_Tomog.mat" "$(sed -n 3p "${CHUNKS_FILE}" | cut -f3)"
    assert_file_lines "${FAIL_FILE}" 0
    assert_file_lines "${DONE_FILE}" 1
    [[ ! -d "${RUN_SCRATCH}/c1" ]]
}
test_collect_failed_chunk_with_nothing_to_retry_is_failed() {
    load_collect_fns
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    touch "${RUN_SCRATCH}/c1/FAIL" "${RUN_SCRATCH}/c1/data/field_retrieval/Tomogram_Field_s1_Tomog.mat"
    collect_finished
    assert_file_lines "${FAIL_FILE}" 1
    assert_file_lines "${CHUNKS_FILE}" 1
}
test_collect_failed_chunk_keeps_scratch_if_copy_fails() {
    load_collect_fns
    rsync() { return 23; }
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    touch "${RUN_SCRATCH}/c1/FAIL" "${RUN_SCRATCH}/c1/data/field_retrieval/a.mat"
    collect_finished
    [[ -d "${RUN_SCRATCH}/c1" ]]
    assert_file_lines "${COPYFAIL_FILE}" 1
}
test_collect_running_job_stays_active() {
    load_collect_fns
    squeue() { echo "100 RUNNING"; }
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    collect_finished
    assert_file_lines "${ACTIVE_FILE}" 1
    [[ ! -f "${RUN_SCRATCH}/c1/JOB_GONE" ]]
}
test_collect_vanished_job_is_failed_only_after_grace_period() {
    load_collect_fns
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    collect_finished                       # first sighting: start the grace timer
    [[ -f "${RUN_SCRATCH}/c1/JOB_GONE" ]]
    [[ ! -f "${RUN_SCRATCH}/c1/FAIL" ]]
    assert_file_lines "${ACTIVE_FILE}" 1
    touch -d '1 minute ago' "${RUN_SCRATCH}/c1/JOB_GONE"
    collect_finished                       # still inside the 2 min grace period
    [[ ! -f "${RUN_SCRATCH}/c1/FAIL" ]]
    touch -d '10 minutes ago' "${RUN_SCRATCH}/c1/JOB_GONE"
    collect_finished                       # grace expired: marked failed
    [[ -f "${RUN_SCRATCH}/c1/FAIL" ]]
    collect_finished                       # next poll processes the failure
    # its only sample never finished -> retried once as a new chunk
    assert_file_lines "${CHUNKS_FILE}" 2
    assert_file_lines "${DONE_FILE}" 1
    assert_file_lines "${FAIL_FILE}" 0
    assert_file_lines "${ACTIVE_FILE}" 0
}
test_collect_late_done_marker_beats_vanished_job() {
    load_collect_fns
    make_active_chunk c1 100 "/x/s1_Tomog.mat"
    collect_finished
    touch "${RUN_SCRATCH}/c1/DONE"         # marker appears late on the shared FS
    collect_finished
    assert_file_lines "${DONE_FILE}" 1
    assert_file_lines "${FAIL_FILE}" 0
}
test_collect_orphaned_chunk_with_recorded_outcome_is_dropped() {
    load_collect_fns
    printf '%s\t%s\n' c1 100 > "${ACTIVE_FILE}"      # scratch dir does not exist
    echo c1 > "${DONE_FILE}"
    collect_finished
    assert_file_lines "${ACTIVE_FILE}" 0
    assert_file_lines "${FAIL_FILE}" 0
}
test_collect_orphaned_chunk_without_outcome_counts_as_failed() {
    load_collect_fns
    printf '%s\t%s\n' c1 100 > "${ACTIVE_FILE}"
    collect_finished
    assert_file_lines "${ACTIVE_FILE}" 0
    assert_file_lines "${FAIL_FILE}" 1
}

run_all "orchestrator.sh"
