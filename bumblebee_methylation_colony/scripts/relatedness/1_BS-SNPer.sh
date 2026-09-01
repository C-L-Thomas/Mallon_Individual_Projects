module load samtools/1.18

#!/bin/bash

BS_SNPPER_DIR="/BS-Snper/BS-Snper"
BAM_DIR="/BS-Snper/sorted_bams"
REF="GCF_910591885.1_iyBomTerr1.2_genomic.fa"
OUT_DIR="/BS-Snper_Outputs/"

mkdir -p "$OUT_DIR"

for BAM in "$BAM_DIR"/*.sorted.bam
do
    FILENAME=$(basename "$BAM")
    SAMPLE=${FILENAME%.sorted.bam}

    echo "Running BS-SNPer for $SAMPLE"

    perl "$BS_SNPPER_DIR/BS-Snper.pl" \
        "$BAM" \
        --fa "$REF" \
        --output "$OUT_DIR/${SAMPLE}.snp.txt" \
        --methcg "$OUT_DIR/${SAMPLE}.meth.cg.txt" \
        --methchg "$OUT_DIR/${SAMPLE}.meth.chg.txt" \
        --methchh "$OUT_DIR/${SAMPLE}.meth.chh.txt" \
        --minhetfreq 0.1 \
        --minhomfreq 0.85 \
        --minquali 15 \
        --mincover 10 \
        --maxcover 1000 \
        --minread2 2 \
        --errorate 0.02 \
        --mapvalue 20 \
        > "$OUT_DIR/${SAMPLE}.vcf" \
        2> "$OUT_DIR/${SAMPLE}.err.log"

done
