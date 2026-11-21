#!/bin/bash
# Check that two BAMs hold the same reads in the same order, ignoring the header (the @PG lines record paths and tool
# versions) and the BGZF compression (which differs between samtools versions).
#
#   bash scripts/compare_bams.sh A.bam B.bam      exit status 0 = identical reads
set -u
[[ $# -eq 2 && -f "$1" && -f "$2" ]] || { echo "usage: $0 A.bam B.bam" >&2; exit 2; }
a=$(samtools view "$1" | md5sum | cut -d' ' -f1)
b=$(samtools view "$2" | md5sum | cut -d' ' -f1)
echo "$a  $1"
echo "$b  $2"
if [[ "$a" == "$b" ]]; then echo "identical reads"; exit 0; fi
echo "reads differ" >&2
exit 1
