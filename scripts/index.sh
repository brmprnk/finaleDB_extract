#!/bin/bash
# Index every BAM in a folder that has no index yet (or an index older than the BAM). scripts/convert.sh already
# indexes each BAM; use this to repair a folder, e.g. after copying BAMs without their .bai files.
#
#   bash scripts/index.sh -d BAMDIR [-j PARALLEL]
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/_common.sh"

BAMDIR= PARALLEL=4
while [[ $# -gt 0 ]]; do
	case "$1" in
		-d) BAMDIR="$2"; shift 2 ;;
		-j) PARALLEL="$2"; shift 2 ;;
		-h|--help) sed -n '2,6p' "$0"; exit 0 ;;
		*) die "unknown argument $1" ;;
	esac
done
[[ -d "${BAMDIR}" ]] || die "folder not found: '${BAMDIR}' (-d)"
command -v samtools > /dev/null || die "samtools not on PATH"

index_bam() {
	local bam="$1"
	if [[ -f "${bam}.bai" && "${bam}.bai" -nt "${bam}" ]]; then return 0; fi
	if samtools index "${bam}"; then echo "indexed ${bam}"; else echo "ERROR: failed to index ${bam}" >&2; return 1; fi
}
export -f index_bam

find "${BAMDIR}" -maxdepth 1 -name '*.bam' -type f | sort | xargs -P "${PARALLEL}" -I {} bash -c 'index_bam "$1"' _ {}
