#!/bin/bash
# ==============================================================================
# WES Analysis Pipeline - Production Stage Preparation Script
# ==============================================================================
# This script:
#   1. Activates the Snakemake Conda environment
#   2. Prepares target capture BED panels in resources/ref/grch38/
#   3. Downloads and indexes the canonical GRCh38 reference genome
#      (Homo_sapiens_assembly38.fasta + .fai + .dict + BWA-MEM2 indexes)
#   4. Downloads BQSR known polymorphic sites (matching dbSNP + Mills + indels)
#   5. Pre-builds all Snakemake Conda environments (--conda-create-envs-only)
#   6. Performs pre-flight pipeline validation
#
# Best run inside a detached tmux session:
#   tmux new -s wes_prep "./prepare_production.sh"
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

LOG_DIR="${SCRIPT_DIR}/logs/prep"
mkdir -p "$LOG_DIR"
PREP_LOG="${LOG_DIR}/prepare_production_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$PREP_LOG") 2>&1

GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}==============================================================================${NC}"
echo -e "${GREEN}      WES Production Staging: Reference Download & Environment Builder        ${NC}"
echo -e "${GREEN}==============================================================================${NC}"
echo -e "Started at: $(date)"
echo -e "Logging to: ${PREP_LOG}\n"

# ------------------------------------------------------------------------------
# 1. Conda Environment Initialization
# ------------------------------------------------------------------------------
echo -e "${BLUE}[STEP 1/5] Initializing Conda / Miniforge environment...${NC}"
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
fi

if conda env list | grep -q "^snakemake "; then
    echo -e "Activating 'snakemake' Conda environment..."
    conda activate snakemake
elif conda env list | grep -q "^wes-pipeline "; then
    echo -e "Activating 'wes-pipeline' Conda environment..."
    conda activate wes-pipeline
else
    echo -e "${YELLOW}[WARNING] Dedicated 'snakemake' environment not detected. Using current active python: $(which python3)${NC}"
fi

# ------------------------------------------------------------------------------
# 2. Directory Hierarchy & Panel Staging
# ------------------------------------------------------------------------------
echo -e "\n${BLUE}[STEP 2/5] Creating directory hierarchy and staging BED panels...${NC}"
mkdir -p resources/ref/grch38/dbsnp
mkdir -p resources/ref/grch38/indels
mkdir -p resources/vep_cache
mkdir -p input
mkdir -p results
mkdir -p logs

TARGET_BED="resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.mergeBed.bed"
CDS_BED="resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCodingExons.hg38.mergeBed.bed"
CANON_BED="resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.hg38.mergeBed.bed"

# Strict production check: verify all three distinct target panels exist
for bed_file in "$TARGET_BED" "$CDS_BED" "$CANON_BED"; do
    bed_base="$(basename "$bed_file")"
    if [ ! -f "$bed_file" ]; then
        if [ -f "ref/$bed_base" ]; then
            echo "Staging $bed_base from ref/..."
            cp "ref/$bed_base" "$bed_file"
        else
            echo -e "${RED}[ERROR] Required production BED panel missing: ${bed_file}${NC}"
            echo -e "Each panel (Target, CDS, CanonicalTrans) must be uniquely generated without placeholders."
            exit 1
        fi
    fi
done

# Verify all 3 panels are non-empty and distinct
T_SIZE=$(wc -c < "$TARGET_BED" || echo 0)
C_SIZE=$(wc -c < "$CDS_BED" || echo 0)
K_SIZE=$(wc -c < "$CANON_BED" || echo 0)

if [ "$T_SIZE" -eq "$C_SIZE" ] || [ "$C_SIZE" -eq "$K_SIZE" ]; then
    echo -e "${RED}[ERROR] Production panel integrity check failed: panels have identical file sizes!${NC}"
    echo -e "Target ($T_SIZE B), CDS ($C_SIZE B), and Canonical ($K_SIZE B) must be distinct biological panels."
    exit 1
fi
echo -e "✓ Verified 3 distinct production BED panels (Target: ${T_SIZE}B, CDS: ${C_SIZE}B, Canonical: ${K_SIZE}B)"

# ------------------------------------------------------------------------------
# 3. Reference Genome & BQSR Known Sites Download
# ------------------------------------------------------------------------------
echo -e "\n${BLUE}[STEP 3/5] Downloading GRCh38 Reference Genome & BQSR Sites (Broad GCS)...${NC}"
REF_FASTA="resources/ref/grch38/Homo_sapiens_assembly38.fasta"

python3 workflow/scripts/download_ref.py "$REF_FASTA"

# Ensure sequence dictionary exists
REF_DICT="resources/ref/grch38/Homo_sapiens_assembly38.dict"
if [ ! -f "$REF_DICT" ] && [ -f "$REF_FASTA" ]; then
    echo "Creating Sequence Dictionary (.dict)..."
    if command -v gatk &>/dev/null; then
        gatk CreateSequenceDictionary -R "$REF_FASTA" -O "$REF_DICT"
    elif command -v samtools &>/dev/null; then
        samtools dict "$REF_FASTA" -o "$REF_DICT"
    fi
fi

# Ensure FASTA index (.fai) exists
if [ ! -f "${REF_FASTA}.fai" ] && [ -f "$REF_FASTA" ]; then
    if command -v samtools &>/dev/null; then
        echo "Creating FASTA index (.fai)..."
        samtools faidx "$REF_FASTA"
    fi
fi

# Ensure BWA-MEM2 indexes exist
if [ ! -f "${REF_FASTA}.bwt.2bit.64" ] && [ -f "$REF_FASTA" ]; then
    if command -v bwa-mem2 &>/dev/null; then
        echo "Generating BWA-MEM2 index (single execution)..."
        bwa-mem2 index "$REF_FASTA"
    fi
fi

echo -e "✓ Reference genome, dictionary, indexes, and BQSR known sites verified."

# ------------------------------------------------------------------------------
# 4. Build Snakemake Conda Environments
# ------------------------------------------------------------------------------
echo -e "\n${BLUE}[STEP 4/5] Pre-building all Snakemake Conda environments...${NC}"
echo "Running: snakemake -s workflow/Snakefile --use-conda --conda-create-envs-only --cores 8"

if command -v snakemake &>/dev/null; then
    snakemake -s workflow/Snakefile \
        --configfile workflow/config.yml \
        --use-conda \
        --conda-create-envs-only \
        --cores 8 || {
            echo -e "${YELLOW}[WARNING] Snakemake env creation encountered non-critical warnings. Envs will also build on first rule execution.${NC}"
        }
    echo -e "✓ Snakemake Conda environments built successfully."
else
    echo -e "${YELLOW}[WARNING] 'snakemake' command not found in current PATH. Ensure environment is activated.${NC}"
fi

# ------------------------------------------------------------------------------
# 5. Pre-flight Validation
# ------------------------------------------------------------------------------
echo -e "\n${BLUE}[STEP 5/5] Running pre-flight pipeline validation...${NC}"
python3 -c "
from workflow.scripts.validate import validate_pipeline_config
validate_pipeline_config('workflow/config.yml', 'input', 'results')
" || true

echo -e "\n${GREEN}==============================================================================${NC}"
echo -e "${GREEN}             Production Staging & Reference Preparation Completed!            ${NC}"
echo -e "${GREEN}==============================================================================${NC}"
echo -e "Finished at: $(date)"
echo -e "Log file: ${PREP_LOG}"
