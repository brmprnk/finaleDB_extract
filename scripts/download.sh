#!/bin/bash
# Download the FinaleDB hg38 fragment files of a list of records.
#
#   bash scripts/download.sh -i IDS -o OUTDIR [-j PARALLEL] [--tbi]
#
#   IDS     plain list of record IDs (EE87786, ...) or a CSV with an ID column (e.g. a sample sheet)
#   OUTDIR  receives <ID>.hg38.frag.tsv.bgz (and .tbi with --tbi)
#
# A file is complete when its size equals the server's Content-Length; complete files are not downloaded again, and
# partial ones are resumed. Records that FinaleDB serves without hg38 fragments (a 28-byte empty BGZF file; their
# metadata API entry is empty too) are not kept: they are listed in OUTDIR/empty_records.txt and the script exits
# with status 2. Failed downloads are listed in OUTDIR/failed_records.txt (exit status 1).
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/_common.sh"

IDS= OUTDIR= PARALLEL=4 TBI=false
while [[ $# -gt 0 ]]; do
	case "$1" in
		-i) IDS="$2"; shift 2 ;;
		-o) OUTDIR="$2"; shift 2 ;;
		-j) PARALLEL="$2"; shift 2 ;;
		--tbi) TBI=true; shift ;;
		-h|--help) sed -n '2,13p' "$0"; exit 0 ;;
		*) die "unknown argument $1" ;;
	esac
done
[[ -f "${IDS}" ]] || die "ID list not found: '${IDS}' (-i)"
[[ -n "${OUTDIR}" ]] || die "no output folder (-o)"
mkdir -p "${OUTDIR}" || die "cannot create ${OUTDIR}"

fetch_one() {  # fetch_one URL DEST -> 0 complete, 1 failed, 2 empty record
	local url="$1" dest="$2" remote local_size
	remote=$(curl -sfIL --retry 5 "${url}" | tr -d '\r' | awk 'tolower($1)=="content-length:" {n=$2} END {print n}')
	[[ -n "${remote}" ]] || { echo "  no response for ${url}" >&2; return 1; }
	if [[ "${remote}" -le ${EMPTY_BGZF_BYTES} ]]; then return 2; fi
	local_size=$([[ -f "${dest}" ]] && wc -c < "${dest}" | tr -d " " || echo 0)
	if [[ "${local_size}" -eq "${remote}" ]]; then echo "  ${dest}: complete, skipped"; return 0; fi
	[[ "${local_size}" -gt "${remote}" ]] && rm -f "${dest}"
	curl -sfL --retry 5 --retry-delay 10 -C - -o "${dest}" "${url}" || { echo "  download failed: ${url}" >&2; return 1; }
	local_size=$([[ -f "${dest}" ]] && wc -c < "${dest}" | tr -d " " || echo 0)
	[[ "${local_size}" -eq "${remote}" ]] || { echo "  ${dest}: ${local_size} of ${remote} bytes" >&2; return 1; }
	echo "  ${dest}: ${remote} bytes"
}

download_record() {
	local id="$1" url dest status
	url=$(finaledb_url "${id}")
	dest="${OUTDIR}/${id}.hg38.frag.tsv.bgz"
	fetch_one "${url}" "${dest}"; status=$?
	if [[ ${status} -eq 2 ]]; then
		echo "  ${id}: FinaleDB has no hg38 fragments for this record (empty file); not kept" >&2
		rm -f "${dest}" "${dest}.tbi"; echo "${id}" >> "${OUTDIR}/empty_records.txt"; return 0
	fi
	if [[ ${status} -eq 0 && ${TBI} == true ]]; then fetch_one "${url}.tbi" "${dest}.tbi" || status=1; fi
	[[ ${status} -eq 0 ]] || echo "${id}" >> "${OUTDIR}/failed_records.txt"
	return 0
}
export -f fetch_one download_record finaledb_url
export OUTDIR TBI FINALEDB_S3 EMPTY_BGZF_BYTES

rm -f "${OUTDIR}/empty_records.txt" "${OUTDIR}/failed_records.txt"
N=$(read_ids "${IDS}" | wc -l)
echo "Downloading ${N} FinaleDB records from ${IDS} to ${OUTDIR} (${PARALLEL} in parallel)"
read_ids "${IDS}" | xargs -P "${PARALLEL}" -I {} bash -c 'download_record "$1"' _ {}

STATUS=0
if [[ -s "${OUTDIR}/failed_records.txt" ]]; then
	echo "$(wc -l < "${OUTDIR}/failed_records.txt") record(s) failed, rerun to resume: ${OUTDIR}/failed_records.txt" >&2; STATUS=1
fi
if [[ -s "${OUTDIR}/empty_records.txt" ]]; then
	sort -o "${OUTDIR}/empty_records.txt" "${OUTDIR}/empty_records.txt"
	echo "$(wc -l < "${OUTDIR}/empty_records.txt") record(s) have no hg38 data on FinaleDB: ${OUTDIR}/empty_records.txt" >&2
	[[ ${STATUS} -eq 0 ]] && STATUS=2
fi
[[ ${STATUS} -eq 0 ]] && echo "All ${N} records complete."
exit ${STATUS}
