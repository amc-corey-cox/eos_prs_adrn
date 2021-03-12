#!/bin/bash
#SBATCH -p c2s4
#SBATCH --job-name=hpc_tmp
#SBATCH --out=hpc_tmp.log
#SBATCH --error=hpc_tmp.err

# These lines generally make bash shell scripts safer
set -euo pipefail        # -e: exit on any error, -u: treat unset variables as error
IFS="`printf '\n\t'`"    # split words only on \n and \t, not space (improves loops)

# Uncomment this for better logging
set -x # Print each command after variable exansion

tmp_dir=$(mktemp -d)
cd "${tmp_dir}"

tools_dir=hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/EOS_PRS/tools

gsutil -q cp "gs://${tools_dir}/bin/tabix" .
chmod +x tabix

PATH=$tmp_dir/:$PATH

mkdir -p gcs_mnt/hdcekabarnesadrn1/
gcsfuse --implicit-dirs hdcekabarnesadrn1 gcs_mnt/hdcekabarnesadrn1/

base_dir=gcs_mnt/hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/b_EA_MEGA/imputed
for i in {22..1}; do
  if [[ -f "${base_dir}/chr${i}.dose.vcf.gz.tbi" ]]; then continue; fi
  echo tabix -f "${base_dir}/chr${i}.dose.vcf.gz"
  tabix -f "${base_dir}/imputed/chr${i}.dose.vcf.gz"
done


