module load bcftools
module load samtools

VCF_DIR="BS-Snper_Outputs/"
BAM_DIR="BS-Snper/sorted_bams/"
OUT_DIR="BCFTools_Formatted/"

REF="GCF_910591885.1_iyBomTerr1.2_genomic.fa"

mkdir -p "$OUT_DIR"

BAM_LIST="$OUT_DIR/bams.txt"
CANDIDATES="$OUT_DIR/transversion_candidates.tsv"
POSITIONS="$OUT_DIR/positions.tsv"
ALLELES="$OUT_DIR/alleles.tsv.gz"

RAW_VCF="$OUT_DIR/joint.raw.vcf.gz"
RENAMED_VCF="$OUT_DIR/joint.renamed.vcf.gz"
FINAL_VCF="$OUT_DIR/bombus_joint.transversions.dp5.vcf.gz"
NEW_NAMES="$OUT_DIR/sample_names.txt"


# ------------------------------------------------------------
# 1. Create the BAM list
# ------------------------------------------------------------

find "$BAM_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.sorted.bam' |
sort > "$BAM_LIST"

echo "BAM files: $(wc -l < "$BAM_LIST")"


# ------------------------------------------------------------
# 2. Extract unique PASS biallelic transversions
# ------------------------------------------------------------
#
# BS-SNPer has already applied:
#   - base quality
#   - mapping quality
#   - minimum and maximum coverage
#   - minimum alternative-read count
#   - allele-frequency thresholds
#
# Here we:
#   - retain BS-SNPer PASS calls
#   - retain biallelic variants
#   - remove transitions
#   - remove duplicate records
#   - remove positions with conflicting alleles

for VCF in "$VCF_DIR"/*.vcf
do

    bcftools view \
        --apply-filters PASS \
        --min-alleles 2 \
        --max-alleles 2 \
        --output-type u \
        "$VCF" |
    bcftools query \
        --format '%CHROM\t%POS\t%REF\t%ALT\n'
done |
awk '
#Remove Bisulphite artifacts
BEGIN {
    FS = OFS = "\t"
}

{
    change = $3 $4

    if (change != "AG" &&
        change != "GA" &&
        change != "CT" &&
        change != "TC") {
        print
    }
}
' |
#Remove duplicate lines
sort -u |
awk '
BEGIN {
    FS = OFS = "\t"
}

{
    position = $1 FS $2
    number[position]++
    record[position] = $0
}

END {
    for (position in number) {
        if (number[position] == 1) {
            print record[position]
        }
    }
}
' |
#Sort in genomic order
sort -k1,1 -k2,2n > "$CANDIDATES"

echo "Candidate transversions: $(wc -l < "$CANDIDATES")"


# ------------------------------------------------------------
# 3. Genotype the candidate sites in all samples
# ------------------------------------------------------------

#Firstly get the coordinates required for examination

awk '
BEGIN {
    OFS = "\t"
}

{
    print $1, $2
}
' "$CANDIDATES" > "$POSITIONS"

awk '
BEGIN {
    OFS = "\t"
}

{
    print $1, $2, $3 "," $4
}
' "$CANDIDATES" |
bgzip -c > "$ALLELES"

tabix \
    --force \
    --sequence 1 \
    --begin 2 \
    --end 2 \
    "$ALLELES"


# Recall SNP at each position

bcftools mpileup \
    --fasta-ref "$REF" \
    --bam-list "$BAM_LIST" \
    --ignore-RG \
    --targets-file "$POSITIONS" \
    --min-MQ 20 \
    --min-BQ 20 \
    --max-depth 1000 \
    --annotate FORMAT/DP,FORMAT/AD \
    --output-type u |
bcftools call \
    --multiallelic-caller \
    --variants-only \
    --constrain alleles \
    --targets-file "$ALLELES" \
    --output-type z \
    --output "$RAW_VCF"


# ------------------------------------------------------------
# 4. Create shorter sample names
# ------------------------------------------------------------

sed -E '
s#.*/##
s/\.deduplicated\.sorted\.bam$//
s/\.sorted\.bam$//
s/_bismark.*//
' "$BAM_LIST" > "$NEW_NAMES"

bcftools reheader \
    --samples "$NEW_NAMES" \
    --output "$RENAMED_VCF" \
    "$RAW_VCF"


# ------------------------------------------------------------
# 5. Set genotypes with depth below 5 to missing
# ------------------------------------------------------------
#Refilter for coverage after bam added values again

bcftools +setGT \
    "$RENAMED_VCF" \
    --output-type z \
    --output "$FINAL_VCF" \
    -- \
    --target-gt q \
    --new-gt . \
    --include 'FMT/DP<5'


# ------------------------------------------------------------
# 6. Index the final VCF
# ------------------------------------------------------------

bcftools index \
    --tbi \
    --force \
    "$FINAL_VCF"

echo
echo "Finished"
echo "Samples:  $(bcftools query -l "$FINAL_VCF" | wc -l)"
echo "Variants: $(bcftools index -n "$FINAL_VCF")"
echo "Final VCF: $FINAL_VCF"
