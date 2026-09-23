module load plink2

# Input joint VCF
VCF="bombus_joint.transversions.dp5.vcf.gz"

# PLINK output directory
OUT_DIR="plink"

mkdir -p "$OUT_DIR"

PLINK_PREFIX="$OUT_DIR/bombus_transversions"
MISSING_PREFIX="$OUT_DIR/bombus_missingness"

# ------------------------------------------------------------
# 1. Convert the VCF to PLINK 2 format
# ------------------------------------------------------------

plink2 \
    --vcf "$VCF" \
    --allow-extra-chr \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --make-pgen \
    --out "$PLINK_PREFIX"

# ------------------------------------------------------------
# 2. Calculate missingness
# ------------------------------------------------------------
#
# .smiss = missingness for each sample
# .vmiss = missingness for each variant

plink2 \
    --pfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --missing \
    --out "$MISSING_PREFIX"

# ------------------------------------------------------------
# 3. Apply quality-control filters
# ------------------------------------------------------------
#
# --geno 0.10:
#   Remove variants missing in more than 10% of samples.
#
# --maf 0.05:
#   Remove variants with a minor allele frequency below 5%.

QC_PREFIX="$OUT_DIR/bombus_qc"

plink2 \
    --pfile "$PLINK_PREFIX" \
    --allow-extra-chr \
    --geno 0.10 \
    --maf 0.05 \
    --make-pgen \
    --out "$QC_PREFIX"


# ------------------------------------------------------------
# 4. Physically thin the SNPs
# ------------------------------------------------------------
#
# PLINK does not estimate LD with fewer than 50 samples.
# Instead, retain SNPs separated by at least 10,000 bp as outlined by b impatiens 
#
# This reduces the effect of large numbers of closely positioned,
# correlated SNPs without estimating LD from the 36 samples.

PRUNED_PREFIX="$OUT_DIR/bombus_pruned"

plink2 \
    --pfile "$QC_PREFIX" \
    --allow-extra-chr \
    --bp-space 10000 \
    --make-pgen \
    --out "$PRUNED_PREFIX"

# ------------------------------------------------------------
# 5. Calculate pairwise KING kinship
# ------------------------------------------------------------
#
# This produces one row for every pair of samples.
# The main result is the KINSHIP column.

KINSHIP_PREFIX="$OUT_DIR/bombus_kinship"

plink2 \
    --pfile "$PRUNED_PREFIX" \
    --allow-extra-chr \
    --make-king-table \
    --out "$KINSHIP_PREFIX"


# ------------------------------------------------------------
# 7. Calculate the genomic relationship matrix
# ------------------------------------------------------------
#
#
# Calculate allele frequencies for the pruned SNPs

FREQUENCY_PREFIX="$OUT_DIR/bombus_pruned_frequencies"

plink2 \
    --pfile "$PRUNED_PREFIX" \
    --allow-extra-chr \
    --freq counts \
    --out "$FREQUENCY_PREFIX"


#
# Calculate the genomic relationship matrix
# 
# --read-freq supplies the allele frequencies explicitly.
# This allows --make-rel to run with fewer than 50 samples.
#
# These frequencies are calculated from the same 36 bees, so the
# result should be treated as an exploratory relationship matrix.

RELATIONSHIP_PREFIX="$OUT_DIR/bombus_relationship"

plink2 \
    --pfile "$PRUNED_PREFIX" \
    --allow-extra-chr \
    --read-freq "$FREQUENCY_PREFIX.acount" \
    --make-rel square \
    --out "$RELATIONSHIP_PREFIX"


# ------------------------------------------------------------
# 8. Report the final outputs
# ------------------------------------------------------------

echo
echo "PLINK analysis completed"
echo "Initial dataset:     ${PLINK_PREFIX}.pgen/.pvar/.psam"
echo "Missingness reports: ${MISSING_PREFIX}.smiss and .vmiss"
echo "QC dataset:          ${QC_PREFIX}.pgen/.pvar/.psam"
echo "LD-pruned dataset:   ${PRUNED_PREFIX}.pgen/.pvar/.psam"
echo "Kinship table:       ${KINSHIP_PREFIX}.kin0"
echo "Relationship matrix: ${RELATIONSHIP_PREFIX}.rel"
echo "Matrix sample order: ${RELATIONSHIP_PREFIX}.rel.id"

# ------------------------------------------------------------
# 9. Calculate principal components
# ------------------------------------------------------------
#
# Calculate the first 10 principal components using the
# physically thinned SNP dataset.
#
# The same allele-frequency file is supplied because the dataset
# contains fewer than 50 samples.
#
# The PCA will show whether bees cluster according to colony.

PCA_PREFIX="$OUT_DIR/bombus_pca"

plink2 \
    --pfile "$PRUNED_PREFIX" \
    --allow-extra-chr \
    --read-freq "$FREQUENCY_PREFIX.acount" \
    --pca 10 \
    --out "$PCA_PREFIX"
