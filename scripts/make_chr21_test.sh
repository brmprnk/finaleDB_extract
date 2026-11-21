#!/bin/bash
# Build the chr21 test set end to end: chr21 fragments of each record -> BAM with the modified Fragmentstein -> index,
# and write record checksums so that a rebuild can be compared with the published one.
#
#   bash scripts/make_chr21_test.sh -g GENOME -o OUTDIR [-i IDS] [-s SETTINGS] [-f FRAGDIR] [-j PARALLEL] [-t THREADS]
#
#   IDS       records (default lists/chr21_test.csv: 10 healthy + 10 BRCA/CRC from the Cristiano cohort)
#   SETTINGS  per-record MAPQ filter, as for scripts/convert.sh (default lists/cris_conversion_settings.csv)
#   FRAGDIR   optional folder with full fragment files from scripts/download.sh; records not found there are read
#             from FinaleDB directly with tabix, which transfers only the chr21 blocks (about 3.5 MB per record)
#   OUTDIR    receives fragments/<ID>.hg38.frag.tsv.bgz (chr21 only), cris/<ID>.hg38.frag.tsv.bam(.bai),
#             checksums.md5 (md5 of the reads of each BAM, header excluded) and sizes.tsv.
#             OUTDIR has the layout of PREPROCESSING_DIR.
#
# If tests/chr21_expected.md5 exists, the checksums are compared with it (exit status 1 on a difference).
# The read names (F1, F2, ...) number the fragments of the chr21 file, so they differ from the names of the same
# reads in a BAM made from the full fragment file; positions, flags, MAPQ, sequences and qualities are the same.
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${SCRIPT_DIR}/.." && pwd)"
. "${SCRIPT_DIR}/_common.sh"

IDS="${REPO}/lists/chr21_test.csv" SETTINGS="${REPO}/lists/cris_conversion_settings.csv" FRAGDIR= OUTDIR= GENOME= PARALLEL=4 THREADS=2
REGION=chr21
while [[ $# -gt 0 ]]; do
	case "$1" in
		-i) IDS="$2"; shift 2 ;;
		-s) SETTINGS="$2"; shift 2 ;;
		-f) FRAGDIR="$2"; shift 2 ;;
		-o) OUTDIR="$2"; shift 2 ;;
		-g) GENOME="$2"; shift 2 ;;
		-j) PARALLEL="$2"; shift 2 ;;
		-t) THREADS="$2"; shift 2 ;;
		-h|--help) sed -n '2,19p' "$0"; exit 0 ;;
		*) die "unknown argument $1" ;;
	esac
done
[[ -f "${IDS}" ]] || die "ID list not found: '${IDS}'"
[[ -n "${OUTDIR}" ]] || die "no output folder (-o)"
[[ -f "${GENOME}" && -f "${GENOME}.fai" ]] || die "reference '${GENOME}' or '${GENOME}.fai' not found (-g)"
for tool in tabix bgzip samtools md5sum; do command -v "${tool}" > /dev/null || die "${tool} not on PATH"; done
mkdir -p "${OUTDIR}/fragments" "${OUTDIR}/cris" || die "cannot create ${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"

# 1. chr21 fragments
subset_record() {
	local id="$1" dest="${OUTDIR}/fragments/$1.hg38.frag.tsv.bgz" local_file="${FRAGDIR:+${FRAGDIR}/$1.hg38.frag.tsv.bgz}"
	if [[ -s "${dest}" && -f "${dest}.ok" ]]; then echo "${id}: chr21 fragments present"; return 0; fi
	rm -f "${dest}" "${dest}.ok"
	if [[ -n "${local_file}" && -f "${local_file}" ]]; then
		gunzip -c "${local_file}" | awk -v c="${REGION}" '$1 == c' | bgzip -c > "${dest}" || { echo "${id}: subsetting failed" >&2; return 1; }
	else
		# -D: read the remote index without saving it; run in a scratch folder in case this htslib saves it anyway
		( cd "${OUTDIR}/fragments" && tabix -D "$(finaledb_url "${id}")" "${REGION}" ) | bgzip -c > "${dest}" \
			|| { echo "${id}: remote tabix failed" >&2; return 1; }
	fi
	if [[ -z "$(gunzip -c "${dest}" | head -c 1)" ]]; then echo "${id}: no ${REGION} fragments" >&2; rm -f "${dest}"; return 1; fi
	touch "${dest}.ok"
	echo "${id}: $(gunzip -c "${dest}" | wc -l) ${REGION} fragments"
}
export -f subset_record finaledb_url
export OUTDIR FRAGDIR REGION FINALEDB_S3

echo "1/3 Extracting ${REGION} fragments of $(read_ids "${IDS}" | wc -l) records"
read_ids "${IDS}" | xargs -P "${PARALLEL}" -I {} bash -c 'subset_record "$1"' _ {} || die "fragment extraction failed for at least one record"

# 2. BAMs
echo "2/3 Converting to BAM"
bash "${SCRIPT_DIR}/convert.sh" -i "${IDS}" -f "${OUTDIR}/fragments" -o "${OUTDIR}/cris" -g "${GENOME}" -s "${SETTINGS}" \
	-j "${PARALLEL}" -t "${THREADS}" --logs "${OUTDIR}/logs" || die "conversion failed"

# 3. checksums and sizes
echo "3/3 Checksums"
: > "${OUTDIR}/checksums.md5"
printf 'ID\tfragments\treads\tfragment_file_bytes\tbam_bytes\n' > "${OUTDIR}/sizes.tsv"
for id in $(read_ids "${IDS}"); do
	bam="cris/${id}.hg38.frag.tsv.bam"
	echo "$(samtools view "${OUTDIR}/${bam}" | md5sum | cut -d' ' -f1)  ${bam}" >> "${OUTDIR}/checksums.md5"
	printf '%s\t%s\t%s\t%s\t%s\n' "${id}" "$(gunzip -c "${OUTDIR}/fragments/${id}.hg38.frag.tsv.bgz" | wc -l)" \
		"$(samtools idxstats "${OUTDIR}/${bam}" | awk '{n+=$3+$4} END {print n}')" \
		"$(wc -c < "${OUTDIR}/fragments/${id}.hg38.frag.tsv.bgz" | tr -d ' ')" "$(wc -c < "${OUTDIR}/${bam}" | tr -d ' ')" >> "${OUTDIR}/sizes.tsv"
done
awk -F'\t' 'NR > 1 {f += $4; b += $5} END {printf "Total: %.1f MB fragment files, %.1f MB BAMs\n", f / 1e6, b / 1e6}' "${OUTDIR}/sizes.tsv"

EXPECTED="${REPO}/tests/chr21_expected.md5"
if [[ -f "${EXPECTED}" ]]; then
	if diff <(sort -k2 "${EXPECTED}") <(sort -k2 "${OUTDIR}/checksums.md5"); then
		echo "Reads identical to ${EXPECTED}"
	else
		echo "ERROR: reads differ from ${EXPECTED}" >&2; exit 1
	fi
else
	echo "No ${EXPECTED} yet: copy ${OUTDIR}/checksums.md5 there to record this build as the reference"
fi
