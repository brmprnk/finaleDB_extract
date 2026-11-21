# Shared helpers for the scripts in scripts/. Sourced, not run.

# FinaleDB hg38 fragment files: <FINALEDB_S3>/<ID>/hg38/<ID>.hg38.frag.tsv.bgz (+ .tbi)
FINALEDB_S3="${FINALEDB_S3:-https://s3.us-east-2.amazonaws.com/finaledb.epifluidlab.cchmc.org/entries}"
# A record without hg38 fragments is served as a bare BGZF end-of-file block of this many bytes
EMPTY_BGZF_BYTES=28

finaledb_url() {  # finaledb_url ID -> URL of the hg38 fragment file
	echo "${FINALEDB_S3}/$1/hg38/$1.hg38.frag.tsv.bgz"
}

read_ids() {  # read_ids FILE -> one record ID per line
	# FILE is a plain list (one ID per line) or a CSV whose header has a column named ID, such as a sample
	# sheet. Lines starting with # and blank lines are skipped; duplicates are dropped.
	awk -F',' '
		/^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
		!seen_header { seen_header = 1; for (i = 1; i <= NF; i++) { h = $i; gsub(/[[:space:]"\r]/, "", h); if (h == "ID") col = i }
		               if (col) next; col = 1 }
		{ id = $col; gsub(/[[:space:]"\r]/, "", id); if (id != "" && !(id in done)) { done[id] = 1; print id } }
	' "$1"
}

die() { echo "ERROR: $*" >&2; exit 1; }
