#!/bin/bash
#SBATCH -p c2s4
#SBATCH --job-name=eos_prs_ext
#SBATCH --out=eos_prs_ext.log
#SBATCH --error=eos_prs_ext.err

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
gsutil -q cp "gs://${data_gs}/tools/bin/tabix" .

# Make tools executable and add to PATH
chmod +x tabix
PATH=$tmp_dir/:$PATH

# Mount Google Storage bucket to local directory
mkdir -p gcs_mnt/hdcekabarnesadrn1/
gcsfuse --implicit-dirs hdcekabarnesadrn1 gcs_mnt/hdcekabarnesadrn1/

# Set input file and location of data files
snp_file=beta_scoring_Trans_EOS_hg38.txt
vcf_dir=gcs_mnt/hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/b_AA/imputed

# Create header file ( disable pipefail ) and extract snps
$( set +o pipefail; zcat "${vcf_dir}/chr1.dose.vcf.gz" | head -n 2 > vcf_header.tmp )
rm -f eos_prs_snps_hg38_extracted.vcf
for i in {1..22}; do
  $( set +o pipefail; zcat "${vcf_dir}/chr${i}.dose.vcf.gz" | grep -m 1 "##contig" >> vcf_header.tmp )
  awk -F '[:_\t]' -v chr="^$i" '$0 ~ chr { print "chr" $1 "\t" $2 "\t" $2 }' "${snp_file}" > extract.tmp
  tabix -R extract.tmp "${vcf_dir}/chr${i}.dose.vcf.gz" >> vcf_snps.tmp
done
$( set +o pipefail; zcat "${vcf_dir}/chr1.dose.vcf.gz" | head -n 20 |\
     grep "#" | grep -v -e "##file" -e "##contig" >> vcf_header.tmp )

# Join vcf header and snps into output vcf file
cat vcf_header.tmp vcf_snps.tmp > eos_prs_snps_hg38_extracted.vcf

# Create name file for working with the snps (I"m not sure we need to do this...)
grep -v "#" eos_prs_snps_hg38_extracted.vcf |\
  awk '{ if ($4 < $5) print $3 "\t" $1 ":" $2 "_" $4 "_" $5; else print $3 "\t" $1 ":" $2 "_" $5 "_" $4; }' \
   > eos_prs_snps_hg38_rename.txt

# Copy files to Google Bucket data directory to keep them
gsutil -q cp eos_prs_snps_hg38_extracted.vcf eos_prs_snps_hg38_rename.txt "gs://${data_gs}/"

# Move log files to hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS/log
