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

tmp_dir=$(mktemp -d)
cd "${tmp_dir}"


data_gs=hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS
gsutil -q cp "gs://${data_gs}/eos_prs_snps_hg38.txt" .
gsutil -q cp "gs://${data_gs}/tools/bin/tabix" .
chmod +x tabix

PATH=$tmp_dir/:$PATH


mkdir -p gcs_mnt/hdcekabarnesadrn1/
gcsfuse --implicit-dirs hdcekabarnesadrn1 gcs_mnt/hdcekabarnesadrn1/

rm -f eos_prs_snps_hg38_extracted.vcf
for i in {1..22}; do
  awk -F '[:_\t]' -v chr="^$i" '$0 ~ chr { print "chr" $1 "\t" $2 "\t" $2 }' \
    eos_prs_snps_hg38.txt > extract.tmp
    tabix -R extract.tmp gcs_mnt/hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/b_AA/imputed/chr$i.dose.vcf.gz >> eos_prs_snps_hg38_extracted.vcf
done

gsutil -q cp eos_prs_snps_hg38_extracted.vcf "gs://${data_gs}/"
