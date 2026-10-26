#!/bin/bash
# test_submit_chunk.sh — tests for staging and job submission in
# scripts/orchestrator.sh (submit_chunk). rsync and sbatch are mocked; the real
# scripts/bash_template.sh is filled in and inspected.
#
#   bash tests/test_submit_chunk.sh

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
SRC_FILE="${REPO_ROOT}/scripts/orchestrator.sh"

# Common setup: a source dir with a background and N samples, mocks, env vars
prepare() {  # [source dir name]
    load_fns sample_index submit_chunk
    SRC="${TOMO_DATA_MOUNT}/proj/${1:-exp}"
    mkdir -p "${SRC}"
    touch "${SRC}/bg001_Tomog.mat"
    TOMO_REPO_DIR="${REPO_ROOT}"
    TOMO_EMAIL="me@example.org"; TOMO_SLURM_PARTITION="gpu"; TOMO_SLURM_GRES="gpu:1"
    TOMO_SLURM_CPUS=8; TOMO_SLURM_MEM=64G; TOMO_SLURM_TIME="24:00:00"
    TOMO_NM="1.337"; TOMO_VERSION="default"
    SBATCH_LOG="${TEST_DIR}/sbatch.log"; : > "${SBATCH_LOG}"
    # rsync -a SRC DEST/ : copy, but fail for any file with FAILME in its name
    rsync() {
        local src="${*: -2:1}" dst="${*: -1}"
        [[ "${src}" == *FAILME* ]] && return 23
        cp "${src}" "${dst}"
    }
    sbatch() {
        echo "$*" >> "${SBATCH_LOG}"
        case "$*" in *--dependency*) echo 9002 ;; *) echo 9001 ;; esac
    }
}
make_samples() {  # names... ; also records them as pending files in sorted order
    local n
    for n in "$@"; do touch "${SRC}/${n}"; done
    ls "${SRC}"/sample*_Tomog.mat | sort > "${STATE_DIR}/pending_files.txt" 2>/dev/null || true
}
csv() {  # names... -> comma-separated full paths
    local out=() n
    for n in "$@"; do out+=("${SRC}/${n}"); done
    (IFS=','; echo "${out[*]}")
}

test_stages_background_and_samples_and_submits() {
    prepare
    make_samples sample001_Tomog.mat sample002_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat sample002_Tomog.mat)"
    local d="${RUN_SCRATCH}/c1/data"
    [[ -f "${d}/bg001_Tomog.mat" && -f "${d}/sample001_Tomog.mat" && -f "${d}/sample002_Tomog.mat" ]]
    [[ -d "${d}/field_retrieval" ]]
    assert_contains "${ACTIVE_FILE}" "^c1	9001$"
    assert_contains "${SUBMIT_TIMES_FILE}" "^c1	[0-9]+$"
}
test_all_backgrounds_are_staged() {
    prepare
    touch "${SRC}/bg002_Tomog.mat"
    make_samples sample001_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    [[ -f "${RUN_SCRATCH}/c1/data/bg001_Tomog.mat" && -f "${RUN_SCRATCH}/c1/data/bg002_Tomog.mat" ]]
}
test_job_script_has_no_unfilled_placeholders() {
    prepare
    make_samples sample001_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    # (the header comment names the literal word __PLACEHOLDERS__; ignore that)
    grep -v '__PLACEHOLDERS__' "${JOBS_DIR}/c1.sh" > "${TEST_DIR}/job_filtered.sh"
    assert_not_contains "${TEST_DIR}/job_filtered.sh" "__[A-Z_0-9]+__"
    bash -n "${JOBS_DIR}/c1.sh"
    assert_contains "${JOBS_DIR}/c1.sh" "DATA_DIR=\"${RUN_SCRATCH}/c1/data\""
    assert_contains "${JOBS_DIR}/c1.sh" "CHUNK_DONE=\"${RUN_SCRATCH}/c1/DONE\""
    assert_contains "${JOBS_DIR}/c1.sh" "CHUNK_FAIL=\"${RUN_SCRATCH}/c1/FAIL\""
    assert_contains "${JOBS_DIR}/c1.sh" "#SBATCH --partition=gpu"
    assert_contains "${JOBS_DIR}/c1.sh" "#SBATCH --time=24:00:00"
}
test_nm_and_version_reach_the_job_script() {
    prepare
    TOMO_NM="1.340"; TOMO_VERSION="alice"
    make_samples sample001_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    assert_contains "${JOBS_DIR}/c1.sh" 'NM="1.340"'
    assert_contains "${JOBS_DIR}/c1.sh" 'VERSION="alice"'
}
test_tiff_index_uses_position_in_whole_source_dir() {
    prepare
    make_samples sample001_Tomog.mat sample002_Tomog.mat sample003_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample002_Tomog.mat sample003_Tomog.mat)"
    local idx="${RUN_SCRATCH}/c1/data/tiff_index.txt"
    assert_file_lines "${idx}" 2
    assert_contains "${idx}" "^sample002_Tomog	2$"
    assert_contains "${idx}" "^sample003_Tomog	3$"
}
test_tiff_offset_counts_samples_in_earlier_chunks_of_same_dir() {
    prepare
    make_samples sample001_Tomog.mat sample002_Tomog.mat sample003_Tomog.mat
    printf 'c0\t%s\t%s\n' "${SRC}" "$(csv sample001_Tomog.mat sample002_Tomog.mat)" >> "${CHUNKS_FILE}"
    printf 'c1\t%s\t%s\n' "${SRC}" "$(csv sample003_Tomog.mat)" >> "${CHUNKS_FILE}"
    submit_chunk c1 "${SRC}" "$(csv sample003_Tomog.mat)"
    assert_contains "${JOBS_DIR}/c1.sh" 'TIFF_OFFSET="2"'
}
test_unstageable_sample_is_skipped_not_fatal() {
    prepare
    make_samples sample001_Tomog.mat sample002_FAILME_Tomog.mat sample003_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat sample002_FAILME_Tomog.mat sample003_Tomog.mat)"
    local d="${RUN_SCRATCH}/c1/data"
    [[ -f "${d}/sample001_Tomog.mat" && -f "${d}/sample003_Tomog.mat" ]]
    [[ ! -f "${d}/sample002_FAILME_Tomog.mat" ]]
    assert_file_lines "${SKIPPED_FILE}" 1
    assert_contains "${SKIPPED_FILE}" "^c1	sample002_FAILME_Tomog.mat	"
    assert_not_contains "${d}/tiff_index.txt" "FAILME"
    assert_file_lines "${ACTIVE_FILE}" 1   # the chunk was still submitted
}
test_fails_without_submitting_when_no_sample_can_be_staged() {
    prepare
    make_samples sample001_FAILME_Tomog.mat
    ! submit_chunk c1 "${SRC}" "$(csv sample001_FAILME_Tomog.mat)"
    assert_file_lines "${ACTIVE_FILE}" 0
    assert_file_lines "${SBATCH_LOG}" 0
}
test_fails_without_a_background() {
    prepare
    rm "${SRC}/bg001_Tomog.mat"
    make_samples sample001_Tomog.mat
    ! submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    assert_file_lines "${ACTIVE_FILE}" 0
}
test_fails_when_background_cannot_be_copied() {
    prepare
    mv "${SRC}/bg001_Tomog.mat" "${SRC}/bg001_FAILME_Tomog.mat"
    make_samples sample001_Tomog.mat
    ! submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    assert_file_lines "${ACTIVE_FILE}" 0
}
test_source_dir_name_with_quote_and_space_yields_valid_script() {
    prepare "O'Brien exp"
    make_samples sample001_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    bash -n "${JOBS_DIR}/c1.sh"
    # sourcing just the assignment must give back the original name
    local line name
    line="$(grep '^BATCH_NAME=' "${JOBS_DIR}/c1.sh")"
    eval "${line}"; name="${BATCH_NAME}"
    assert_eq "O'Brien exp" "${name}"
}
test_start_notification_is_submitted_only_once() {
    prepare
    make_samples sample001_Tomog.mat sample002_Tomog.mat
    submit_chunk c1 "${SRC}" "$(csv sample001_Tomog.mat)"
    submit_chunk c2 "${SRC}" "$(csv sample002_Tomog.mat)"
    assert_eq 1 "$(grep -c -- '--dependency' "${SBATCH_LOG}")"
    assert_eq 9002 "$(cat "${STATE_DIR}/start_notify_jid")"
    assert_file_lines "${ACTIVE_FILE}" 2
}

run_all "orchestrator.sh submit_chunk"
