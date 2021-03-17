#!/bin/bash
#SBATCH -p c2s4
#SBATCH --job-name=eos_plink
#SBATCH --out=eos_plink.log
#SBATCH --error=eos_plink.err

# These lines generally make bash shell scripts safer
set -euo pipefail        # -e: exit on any error, -u: treat unset variables as error
IFS="`printf '\n\t'`"    # split words only on \n and \t, not space (improves loops)

# Uncomment this for better logging
set -x # Print each command after variable exansion

# Create tmp directory to work in
tmp_dir=$(mktemp -d)
cd "${tmp_dir}"


# Copy data and tools to working directory
data_gs=hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS
gsutil -q cp "gs://${data_gs}/beta_scoring_Trans_EOS_hg38.txt" .
gsutil -q cp "gs://${data_gs}/tools/bin/plink" .

# Make tools executable and add to PATH
chmod +x plink
PATH=$tmp_dir/:$PATH

# Mount Google Storage bucket to local directory
mkdir -p gcs_mnt/hdcekabarnesadrn1/
gcsfuse --implicit-dirs hdcekabarnesadrn1 gcs_mnt/hdcekabarnesadrn1/

# Set input file and location of data files
data_dir=gcs_mnt/hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS

# Run plink to get Polygenic Trait Scores (PTS)
plink --vcf "${data_dir}/eos_prs_snps_hg38_extracted.vcf" \
  --update-name "${data_dir}/eos_prs_snps_hg38_rename.txt" \
  --double-id \
  --score beta_scoring_Trans_EOS_hg38.txt \
  --out beta_score_eos_prs.txt

# Copy files to Google Bucket data directory to keep them
gsutil -q cp beta_score_eos_prs.txt "gs://${data_gs}/"

# Move log files to hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS/log
