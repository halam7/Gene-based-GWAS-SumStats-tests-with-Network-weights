#!/bin/bash
#SBATCH --job-name=sats_gene
#SBATCH --array=1-127
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sats_gene_%A_%a.out
#SBATCH --error=logs/sats_gene_%A_%a.err

cd /home/halam7/epid_meidaw_chi_link/snp_network

# path for plink12
export PATH="/home/halam7:$PATH"


# Load R module (adjust if needed for your cluster)
module load R/4.2.2-foss-2022b

# Print the R version being used
echo "R version check:"
Rscript --version

# Set the trait file manually or pass it in with sbatch --export
TRAIT_FILE=${TRAIT_FILE:-"data/processed/sumstats_weights2/traitX_weight.txt"}

echo "Running chunk $SLURM_ARRAY_TASK_ID for trait file $TRAIT_FILE"

# Run the chunk
Rscript scripts/17_simulation_genes_parallel.R "$TRAIT_FILE" "$SLURM_ARRAY_TASK_ID"
