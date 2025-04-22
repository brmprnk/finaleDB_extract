#!/bin/bash

VERSION=2023.4a

INPUT_DIR="jiang_frag"
OUTPUT_DIR="jiang"

# Accept number of processes as a command-line argument or default to 1
NUM_PROCESSES=${1:-1}

process_file() {

    CHROM_SIZES="reference_genome/hg38.analysisSet.fa.fai"
    GENOME="reference_genome/hg38.analysisSet.fa"

    READ_LENGTH=101
    MAP_QUALITY_FILTER=5
    MAP_QUALITY_DEFAULT=5
    BASE_QUALITY="F"
    SORT=true
    INCOMPLETE_N=false
    THREADS=4

    INPUT="$1"
    OUTPUT="${OUTPUT_DIR}/$(basename "${INPUT/.bgz/.bam}")"
    TMPBED="${OUTPUT/.bam/.tmp.bed}"

    if [[ -f "${OUTPUT}" ]]; then
        echo "Output file ${OUTPUT} already exists, skipping this step"
        return
    fi

    echo "Processing ${INPUT} -> ${OUTPUT}"
    echo "Unzipping, converting FinaleDB file ${INPUT} and accepting mapping quality>${MAP_QUALITY_FILTER} to BED file format ..."
    echo "Creating paired reads and add fragment length etc..."

    # Script starts
    gunzip -c "$INPUT" | \
    awk -v mqf=${MAP_QUALITY_FILTER} '$4>mqf' | \
    awk -F"\t" -v OFS="\t" '{{k=$4; l=$5; $4="F"++i; $5=k; $6=l}}1' | \
    awk -v r="${READ_LENGTH}" 'BEGIN{FS=OFS="\t"}
        { x= $3; f= $3-$2;
        if ($2+r < $3) $3= $2+r;
        $6= "+";
        $7= f;
        }
        {print}
        { if (x-r > $2) $2= x-r;
        $3= x;
        $6= "-";
        $7= f;
        }
        {print}' > "${TMPBED}"

    TMPSAM=${TMPBED/.bed/.sam}
    echo "Create intermediate SAM file ${TMPSAM} from ${TMPBED} for chromosome_sizes:${CHROM_SIZES} and default map quality: ${MAP_QUALITY_DEFAULT} ..."
    bedToBam -i "${TMPBED}" -mapq "${MAP_QUALITY_DEFAULT}" -g "${CHROM_SIZES}" | \
    samtools view -@ ${THREADS} - | cut -f-8 | \
    paste - <(cut -f5,7 "${TMPBED}") | \
    awk 'BEGIN{FS=OFS="\t"} { print $1,$2,$3,$4,$9,$6,$7,$8,$10 }' > "${TMPSAM}" \

    # extract header, add fragmentstein info
    echo "Creating header for ${TMPSAM} ..."
    HEADER_INFO="@PG	ID:fragmentstein	PN:fragmentstein	VN:${VERSION}	CL:$0 -i ${INPUT} -o ${OUTPUT} -g ${GENOME} -c ${CHROM_SIZES} -r ${READ_LENGTH} -qf ${MAP_QUALITY_FILTER} -qd ${MAP_QUALITY_DEFAULT} -s ${SORT} -t ${THREADS}"
    bedToBam -i "${TMPBED}" -mapq "${MAP_QUALITY_DEFAULT}" -g "${CHROM_SIZES}" | \
    samtools view -H - | \
    awk -v hi="${HEADER_INFO}" 'NR==2 {print hi}1' > "${TMPSAM}.header"

    echo "Get sequences from reference fasta ${GENOME} with read length ${READ_LENGTH} for intermediate files ${TMPBED} and ${TMPSAM} ... "
    echo "if read fwd->flag=99, if read rev -> flag=147, set mate start, set mate chr to same as this read and set base quality to F ... "
    # log command: for debug
    # echo "bedtools getfasta -name -fi \"${GENOME}\" -bed \"${TMPBED}\" | grep -v \">\" | tr '[:lower:]' '[:upper:]' | paste \"${TMPSAM}\" - | awk -v r=\"${READ_LENGTH}\" 'BEGIN{FS=OFS=\"\t\"} {if (\$2==0) {\$2=99; \$8=\$4+\$9-r} else {\$2=147; \$8=\$4-\$9+r; \$9=\$9*-1}} {\$7=\"=\"} {\$11=\$10} {gsub(/./, \"F\", \$11)}{print}' | cat \"${TMPSAM}.header\" - | samtools view -@ ${THREADS} -bo \"${OUTPUT}\""

    bedtools getfasta -name -fi "${GENOME}" -bed "${TMPBED}" | grep -v ">" | tr '[:lower:]' '[:upper:]' | \
    paste "${TMPSAM}" - | \
    awk -v r="${READ_LENGTH}" -v n="${INCOMPLETE_N}" -v bq="${BASE_QUALITY}" 'BEGIN{FS=OFS="\t"} {
            if ($2 == 0) {
                $2= 99;
                if ($4+$9-r > $4) $8= $4+$9-r; else $8= $4;
            } else {
                $2= 147;
                if ($4-$9+r < $4) $8= $4-$9+r; else $8= $4;
                $9= -$9;
            }
            $7= "="; $11= $10;
        } n=="true" {gsub(/B|R|S|W|K|M|Y|H|V|D/, "N", $10)} {gsub(/./, bq, $11)} {print}' | \
    cat "${TMPSAM}.header" - | \
    samtools view -@ ${THREADS} -bo "${OUTPUT}" - \
    && rm "${TMPBED}" && rm "${TMPSAM}" && rm "${TMPSAM}.header"

    if [[ -f "${TMPSAM}" ]]; then
        echo "ERROR: Something went wrong while converting ${TMPBED} and ${TMPSAM} to the final ${OUTPUT} SAM file"
        echo "Investigate the intermediate file with this commands:"
        echo "bedtools getfasta -name -fi \"${GENOME}\" -bed \"${TMPBED}\" | grep -v \">\" | tr '[:lower:]' '[:upper:]' | paste \"${TMPSAM}\" - | awk -v r=\"${READ_LENGTH}\" -v n=\"$INCOMPLETE_N\" 'BEGIN{FS=OFS=\"\\t\"} { if (\$2 == 0) { \$2= 99; if (\$4+\$9-r > \$4) \$8= \$4+\$9-r; else \$8= \$4; } else { \$2= 147; if (\$4-\$9+r < \$4) \$8= \$4-\$9+r; else \$8= \$4; \$9= -\$9; }; \$7= \"=\"; \$11= \$10; } n==\"true\" {gsub(/B|R|S|W|K|M|Y|H|V|D/, \"N\", \$10)} {gsub(/./, \"F\", \$11)} {print}' | cat \"${TMPSAM}.header\" - > ${TMPSAM}2"
        echo "samtools view -@ ${THREADS} -bo ${OUTPUT} ${TMPSAM}2"
        exit 1
    fi

    # Sorting can be done additionally
    if [[ $SORT == true ]]; then
        UNSORTED=${OUTPUT/.bam/.unsorted.bam}
        mv "${OUTPUT}" "${UNSORTED}"
        samtools sort -@ ${THREADS} "${UNSORTED}" -o "${OUTPUT}"
        samtools index -@ ${THREADS} "${OUTPUT}"
        rm "${UNSORTED}"
    fi

    OUTPUT_BW=${OUTPUT/.bam/.bw}
    bamCoverage -b "${OUTPUT}" -o "${OUTPUT_BW}" -p 1

    echo "Finished ${OUTPUT}"
}

export -f process_file
export OUTPUT_DIR MAP_QUALITY_FILTER

find "${INPUT_DIR}" -name "*.bgz" -type f | sort -r | xargs -n 1 -P "${NUM_PROCESSES}" -I {} bash -c 'process_file "$@"' _ {}