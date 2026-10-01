#!/bin/bash
# ==============================================================================
# WES Analysis Pipeline - High-Throughput HPC Slurm Cluster Runner
# 4 Nodes x 80 Cores = 320 Cores Total
# ==============================================================================

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# ------------------------------------------------------------------------------
# 1. Colors & Banner
# ------------------------------------------------------------------------------
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${GREEN}==============================================================================${NC}"
echo -e "${GREEN}        WES Analysis Pipeline - 320-Core HPC Cluster Orchestrator            ${NC}"
echo -e "${GREEN}==============================================================================${NC}"

# ------------------------------------------------------------------------------
# 2. Conda / Miniforge Environment Initialization
# ------------------------------------------------------------------------------
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi

if conda env list | grep -q "^snakemake "; then
    echo -e "${BLUE}[INFO] Activating 'snakemake' Conda environment...${NC}"
    conda activate snakemake
else
    echo -e "${YELLOW}[WARNING] Conda environment 'snakemake' not found.${NC}"
    echo -e "${YELLOW}Create it on the head node using:${NC}"
    echo "  conda create -n snakemake -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm"
fi

# ------------------------------------------------------------------------------
# 3. Cluster Availability & Pre-flight Diagnostics
# ------------------------------------------------------------------------------
echo -e "\n${BLUE}[INFO] Checking cluster status...${NC}"
if command -v sinfo &>/dev/null; then
    sinfo -o "%10P %10N %10c %20C %15m %10e %10T"
    echo -e "\n${BLUE}[INFO] Available Cores (Allocated/Idle/Other/Total):${NC}"
    sinfo -o "%C"
else
    echo -e "${YELLOW}[WARNING] 'sinfo' not found. Ensure Slurm client tools are available.${NC}"
fi

# ------------------------------------------------------------------------------
# 4. Parse Command Line Arguments
# ------------------------------------------------------------------------------
INPUT_DIR=""
OUT_DIR="results"
CONFIG_FILE="workflow/config.yml"
JOBS=320
DEFAULT_MEM=4000
PRE_CREATE_ENVS=false
DRY_RUN=false
DEEPVARIANT=false
EXTRA_ARGS=()

show_usage() {
    cat << EOF
Usage: ./run_cluster.sh -i <input_dir> -o <output_dir> [OPTIONS]

Options:
  -i, --inputdir DIR        Path to raw FASTQ directory (inside /home) [Required]
  -o, --outdir DIR          Path to output directory (inside /home) [Default: results]
  -c, --configfile FILE     Path to config YAML [Default: workflow/config.yml]
  -j, --jobs N              Max concurrent Slurm jobs across cluster [Default: 320]
  -m, --default-mem MB      Default memory in MB for fallback jobs [Default: 4000]
  --deepvariant             Enable dual variant calling (GATK + DeepVariant)
  --pre-create-envs         Pre-build all Conda envs on head node (Compute nodes have no internet)
  -n, --dry-run             Perform dry run DAG resolution without job submission
  -h, --help                Show this message and exit

Examples:
  # Standard 320-core run:
  ./run_cluster.sh -i input/cohort_01 -o results/cohort_01

  # Run inside a resilient tmux session:
  tmux new -s wes_run
  ./run_cluster.sh -i input/cohort_01 -o results/cohort_01
  # (Detach with Ctrl+B then D; reattach later with: tmux attach -t wes_run)
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--inputdir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--outdir)
            OUT_DIR="$2"
            shift 2
            ;;
        -c|--configfile)
            CONFIG_FILE="$2"
            shift 2
            ;;
        -j|--jobs)
            JOBS="$2"
            shift 2
            ;;
        -m|--default-mem)
            DEFAULT_MEM="$2"
            shift 2
            ;;
        --deepvariant)
            DEEPVARIANT=true
            shift
            ;;
        --pre-create-envs)
            PRE_CREATE_ENVS=true
            shift
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

if [ -z "$INPUT_DIR" ]; then
    echo -e "${RED}[ERROR] Input directory (-i / --inputdir) is required.${NC}"
    show_usage
    exit 1
fi

# Ensure paths are within /home (required for compute node visibility)
REAL_INPUT="$(realpath "$INPUT_DIR" 2>/dev/null || echo "$INPUT_DIR")"
REAL_OUT="$(realpath -m "$OUT_DIR" 2>/dev/null || echo "$OUT_DIR")"

if [[ "$REAL_INPUT" != /home/* ]] || [[ "$REAL_OUT" != /home/* ]]; then
    echo -e "${YELLOW}[WARNING] Compute nodes can ONLY access /home.${NC}"
    echo -e "${YELLOW}Input: $REAL_INPUT | Output: $REAL_OUT${NC}"
    echo -e "${YELLOW}Ensure all input/output files reside within /home before running on compute nodes.${NC}"
fi

# ------------------------------------------------------------------------------
# 5. Pre-create Conda Environments (Head Node has internet; Compute nodes do not)
# ------------------------------------------------------------------------------
if [ "$PRE_CREATE_ENVS" = true ]; then
    echo -e "\n${BLUE}[INFO] Pre-building Conda environments on head node...${NC}"
    snakemake -s workflow/Snakefile --configfile "$CONFIG_FILE" --sdm conda --conda-create-envs-only
    echo -e "${GREEN}[SUCCESS] All Conda environments pre-built and cached under .snakemake/conda!${NC}\n"
fi

# ------------------------------------------------------------------------------
# 6. Execute Snakemake on Slurm Cluster
# ------------------------------------------------------------------------------
CMD=(
    python3 wes_pipeline.py run
    -i "$INPUT_DIR"
    -o "$OUT_DIR"
    -c "$CONFIG_FILE"
    --slurm
    --jobs "$JOBS"
    --default-mem "$DEFAULT_MEM"
)

if [ "$DEEPVARIANT" = true ]; then
    CMD+=(--deepvariant)
fi

if [ "$DRY_RUN" = true ]; then
    CMD+=(-n)
fi

if [ ${#EXTRA_ARGS[@]} -gt 0 ]; then
    CMD+=("${EXTRA_ARGS[@]}")
fi

echo -e "\n${GREEN}[INFO] Launching WES pipeline across cluster (Jobs: $JOBS)...${NC}"
echo -e "${BLUE}Command: ${CMD[*]}${NC}\n"

"${CMD[@]}"

echo -e "\n${GREEN}==============================================================================${NC}"
echo -e "${GREEN}                Pipeline run finished successfully!                           ${NC}"
echo -e "${GREEN}==============================================================================${NC}"
