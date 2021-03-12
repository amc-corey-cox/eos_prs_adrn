#!/bin/bash
# These lines generally make bash shell scripts safer
set -euo pipefail        # -e: exit on any error, -u: treat unset variables as error
IFS="`printf '\n\t'`"    # split words only on \n and \t, not space (improves loops)

# Uncomment this for better logging
set -x # Print each command after variable exansion

for i in {1..22}; do
  awk -F '[:_\t]' -v chr="^$i" '$0 ~ chr { print $1 ":" $2 "-" $2 }' \
    BCX_trans_ethnic_risk_scores/beta_scoring_Trans_EOS.bed.lifted.hg38 \
    > extract.tmp
    tabix -T extract.tmp /gcs_mnt/hdcekabarnesadrn1/ADRN_MEGA_META_SEVERITY/b_AA/imputed/chr$i.dose.vcf.gz > tmp.vcf
done
  
