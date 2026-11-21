#!/bin/bash
# Convert downloaded FinaleDB fragment files to sorted, indexed BAMs with the modified Fragmentstein.
#
#   bash scripts/convert.sh -i IDS -f FRAGDIR -o OUTDIR -g GENOME [-s SETTINGS] [-j PARALLEL] [-t THREADS]
#                           [--temp DIR] [--logs DIR] [--bigwig] [--force]
#
#   IDS       plain list of record IDs or a CSV with an ID column
#   FRAGDIR   holds <ID>.hg38.frag.tsv.bgz (scripts/download.sh)
#   OUTDIR    receives <ID>.hg38.frag.tsv.bam and .bam.bai
#   GENOME    hg38.analysisSet.fa with its .fai next to it
#   SETTINGS  optional CSV ID,map_quality_filter with the per-record MAPQ filter (records not listed: 5); pass
#             lists/cris_conversion_settings.csv to reproduce the original Cristiano BAMs
#   PARALLEL  records converted at the same time (default 1); THREADS samtools threads per record (default 4)
#   --temp    folder for the temporary BED/SAM files (roughly 12 GB per 40 million fragments, an estimate); default OUTDIR
#   --logs    folder for one log per record (default logs/convert)
#
# A record is done when its BAM and index exist and pass samtools quickcheck; done records are skipped unless --force.
# Failed records are listed in <logs>/failed_records.txt (exit status 1).
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${SCRIPT_DIR}/_common.sh"
FRAGMENTSTEIN="${SCRIPT_DIR}/../fragmentstein/fragmentstein_finaledb.sh"

IDS= FRAGDIR= OUTDIR= GENOME= SETTINGS= PARALLEL=1 THREADS=4 TEMP= LOGDIR=logs/convert BIGWIG= FORCE=false
while [[ $# -gt 0 ]]; do
	case "$1" in
		-i) IDS="$2"; shift 2 ;;
		-f) FRAGDIR="$2"; shift 2 ;;
		-o) OUTDIR="$2"; shift 2 ;;
		-g) GENOME="$2"; shift 2 ;;
		-s) SETTINGS="$2"; shift 2 ;;
		-j) PARALLEL="$2"; shift 2 ;;
		-t) THREADS="$2"; shift 2 ;;
		--temp) TEMP="$2"; shift 2 ;;
		--logs) LOGDIR="$2"; shift 2 ;;
		--bigwig) BIGWIG=--bigwig; shift ;;
		--force) FORCE=true; shift ;;
		-h|--help) sed -n '2,19p' "$0"; exit 0 ;;
		*) die "unknown argument $1" ;;
	esac
done
[[ -f "${IDS}" ]] || die "ID list not found: '${IDS}' (-i)"
[[ -d "${FRAGDIR}" ]] || die "fragment folder not found: '${FRAGDIR}' (-f)"
[[ -n "${OUTDIR}" ]] || die "no output folder (-o)"
[[ -f "${GENOME}" && -f "${GENOME}.fai" ]] || die "reference '${GENOME}' or its index '${GENOME}.fai' not found (-g)"
[[ -z "${SETTINGS}" || -f "${SETTINGS}" ]] || die "settings file not found: ${SETTINGS}"
for tool in samtools bedtools bedToBam awk gunzip; do command -v "${tool}" > /dev/null || die "${tool} not on PATH"; done
mkdir -p "${OUTDIR}" "${LOGDIR}" || die "cannot create ${OUTDIR} or ${LOGDIR}"
[[ -n "${TEMP}" ]] && { mkdir -p "${TEMP}" || die "cannot create ${TEMP}"; }

convert_record() {
	local id="$1" in out log qf
	in="${FRAGDIR}/${id}.hg38.frag.tsv.bgz"
	out="${OUTDIR}/${id}.hg38.frag.tsv.bam"
	log="${LOGDIR}/${id}.log"
	if [[ ${FORCE} != true && -f "${out}.bai" ]] && samtools quickcheck "${out}" 2> /dev/null; then
		echo "${id}: done, skipped"; return 0
	fi
	if [[ ! -f "${in}" ]]; then
		echo "${id}: ${in} not found" >&2; echo "${id}" >> "${LOGDIR}/failed_records.txt"; return 0
	fi
	qf=5
	if [[ -n "${SETTINGS}" ]]; then
		qf=$(awk -F',' -v id="${id}" '$1 == id {print $2; exit}' "${SETTINGS}" | tr -d '\r'); qf=${qf:-5}
	fi
	rm -f "${out}" "${out}.bai"
	echo "${id}: converting (MAPQ filter ${qf}), log ${log}"
	local start=${SECONDS}
	if bash "${FRAGMENTSTEIN}" -i "${in}" -o "${out}" -g "${GENOME}" -qf "${qf}" -qd "${qf}" -t "${THREADS}" \
		${TEMP:+--temp "${TEMP}"} ${BIGWIG} > "${log}" 2>&1; then
		echo "${id}: finished in $(( (SECONDS - start) / 60 )) min"
	else
		echo "${id}: FAILED, see ${log}" >&2; echo "${id}" >> "${LOGDIR}/failed_records.txt"
	fi
	return 0
}
export -f convert_record
export FRAGDIR OUTDIR GENOME SETTINGS THREADS TEMP LOGDIR BIGWIG FORCE FRAGMENTSTEIN

rm -f "${LOGDIR}/failed_records.txt"
N=$(read_ids "${IDS}" | wc -l)
echo "Converting ${N} records from ${FRAGDIR} to ${OUTDIR} (${PARALLEL} in parallel, ${THREADS} threads each)"
read_ids "${IDS}" | xargs -P "${PARALLEL}" -I {} bash -c 'convert_record "$1"' _ {}

if [[ -s "${LOGDIR}/failed_records.txt" ]]; then
	echo "$(wc -l < "${LOGDIR}/failed_records.txt") record(s) failed: ${LOGDIR}/failed_records.txt" >&2; exit 1
fi
echo "All ${N} records converted."
