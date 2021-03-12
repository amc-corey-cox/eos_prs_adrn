# ====================================== #
#   CHRIS AREHART GUIDE TO A PLINK PRS   #
# ====================================== #

###### QC BEFORE PRS ######
# You will want to do some QC... this article describes QC practices for PRS
# https://www.nature.com/articles/s41596-020-0353-1

# Step 1: R2 > 0.7 or 0.8
# I'll let you do this filtering however you like
# Step 2: Hardy-Weinberg Equilibrium P > 1 × 10−6 using --hwe 0.0000010
# Step 2: Genotyping rate > 0.99 using --geno 0.01
# Step 2: Minor allele frequency (MAF) > 1% using --maf 0.01
plink --bfile myplink \
  --hwe 0.000001 --geno 0.01 --maf 0.01 \
  --make-bed --out hwe_geno_maf
# Step 3: sample missingness and heterozygosity QC
# Step 3: LD prune the data and check for heterozygosity within 3 standard deviations of the mean
plink --bfile hwe_geno_maf --indep 50 5 2 --out hwe_geno_maf_LDpruned_SNPs
# Step 3: Subset the files to the pruned SNPs:
plink --bfile hwe_geno_maf --extract hwe_geno_maf_LDpruned_SNPs.prune.in --make-bed --out hwe_geno_maf_LDpruned
# Step 3: make heterozygosity and missingness files
plink --bfile hwe_geno_maf_LDpruned --het --out hwe_geno_maf_LDpruned
plink --bfile hwe_geno_maf --missing --out hwe_geno_maf
# Step 3: make plots for removal
cat << 'EOF' >c_imiss-vs-het.R
library(data.table)
library("geneplotter")
library("ggplot2")
args <- commandArgs(trailingOnly = TRUE)
inputpathname_het = "./hwe_geno_maf_LDpruned"
inputpathname_imiss = "./hwe_geno_maf"
outputpathname = "imiss-vs-het.pdf"
hetlinesd = 3
missingnessline = 0.02

print(paste("plotting with sd =", hetlinesd, "and miss =s", missingnessline))
imiss = read.table(paste(inputpathname_imiss, ".imiss", sep=""), h=T)
imiss$logF_MISS = log10(imiss[,6]+.00001) # add small value so we dont log10(zero)
het = read.table(paste(inputpathname_het, ".het", sep=""), h=T)
het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.
toPlot <- as.data.frame(cbind(imiss$logF_MISS, het$meanHet))
names(toPlot) <- c("logF_MISS","meanHet")
toPlot$pass_fail <- "PASS"
lowerBound <- mean(het$meanHet) - (hetlinesd*sd(het$meanHet))
upperBound <- mean(het$meanHet) + (hetlinesd*sd(het$meanHet))
toPlot$pass_fail[which(toPlot$logF_MISS >= log10(missingnessline))] <- "FAIL"
summary(as.factor(toPlot$pass_fail))
toPlot$pass_fail[which(toPlot$meanHet < lowerBound)] <- "FAIL"
summary(as.factor(toPlot$pass_fail))
toPlot$pass_fail[which(toPlot$meanHet > upperBound)] <- "FAIL"
toPlot$pass_fail <- as.factor(toPlot$pass_fail)
summary(toPlot$pass_fail)

pdf(outputpathname)
ggplot(toPlot, aes(x=logF_MISS, y=meanHet, color=pass_fail, fill=pass_fail, shape = pass_fail)) +
  geom_point(show.legend = T, alpha = 0.5) + 
  scale_fill_manual(values=rev(c("#D7191C","#2B83BA")))+
  scale_color_manual(values=rev(c("#D7191C","#2B83BA")))+
  theme_light(base_size = 16)+
  xlab("Log10 Proportion of missing genotypes") + ylab("Heterozygosity rate") + labs(title=paste0(""))+
  theme(legend.position="bottom",legend.title = element_blank())+
  geom_vline(aes(xintercept=log10(missingnessline)))+
  geom_hline(aes(yintercept=lowerBound))+
  geom_hline(aes(yintercept=upperBound))+
  scale_shape_manual(values=c(4,16))
ggplot(toPlot, aes(x=meanHet, color=pass_fail, fill=pass_fail, shape = pass_fail)) +
  # geom_point(show.legend = T, alpha = 0.5) +
  geom_histogram(alpha = 0.15, aes(y = ..count..), position = 'identity',show.legend = T, bins = 46) + 
  scale_fill_manual(values=rev(c("#D7191C","#2B83BA")))+
  scale_color_manual(values=rev(c("#D7191C","#2B83BA")))+
  theme_light(base_size = 16)+
  xlab("Heterozygosity rate") + ylab("N") + labs(title=paste0("ALL: Pass=",length(which(toPlot$pass_fail == "PASS")),", Fail=",length(which(toPlot$pass_fail == "FAIL") ) ))+
  theme(legend.position="bottom",legend.title = element_blank())+
  # geom_vline(aes(xintercept=log10(missingnessline)))+
  geom_vline(aes(xintercept=mean(het$meanHet)),linetype="dotted")+
  geom_vline(aes(xintercept=lowerBound))+
  geom_vline(aes(xintercept=upperBound))
dev.off()

FIDlist <- c()
IIDlist <- c()
for(i in 1:nrow(het)){
  if(het$meanHet[i] < mean(het$meanHet) - (hetlinesd*sd(het$meanHet)) | het$meanHet[i] > mean(het$meanHet) + (hetlinesd*sd(het$meanHet))){
  # if(het$meanHet[i] > mean(het$meanHet) + (hetlinesd*sd(het$meanHet))){
    FIDlist <- c(FIDlist,paste0(het$FID[i]))
    IIDlist <- c(IIDlist,paste0(het$IID[i]))
  }
}
cbind(FIDlist,IIDlist)
write.table(cbind(FIDlist,IIDlist),"removeQCfailures.txt",col.names=F,row.names=F,quote=F)
EOF
Rscript c_imiss-vs-het.R

# Step 3: sample missingness < 0.02, 
awk '($5 >= 0.02 && NR>1) {print $2}' hwe_geno_maf.lmiss  > high_miss_SNPs.txt
plink --bfile hwe_geno_maf --exclude high_miss_SNPs.txt --make-bed --out hwe_geno_maf_rm_high_miss_SNPs
# Step 3: remove heterozygosity failures using removeQCfailures.txt
head removeQCfailures.txt
plink --bfile hwe_geno_maf_rm_high_miss_SNPs --remove removeQCfailures.txt --make-bed --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures

# Step 4: remove duplicate SNPs
echo "# *-*-*-* remove duplicate SNPs"
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures --list-duplicate-vars ids-only suppress-first --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures --exclude hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures.dupvar --make-bed --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar

# Step 4: Prioritize some SNPs with a duplicated location (retain the simple biallelic one at the duplicated loci, else remove both)
echo "# *-*-*-* Prioritize some SNPs with the duplicated location (retain the simple biallelic one of the duplicated pair)"
awk '{print $1 ":" $4}' hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar.bim | sort | uniq -d > loc_dup.txt
cat << 'EOF' >rem_loc_dup.R
library(data.table)
loc_dup <- fread("loc_dup.txt",header=F)
colnames(loc_dup) <- c("pos")
if(nrow(loc_dup)>0){
  remDupVar <- read.table("hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar.bim",header=F)
  remDupVar$pos <- paste0(remDupVar$V1,':',remDupVar$V4)
  merged <- merge(loc_dup, remDupVar, by="pos")
  list_to_keep <- c()
  positions_to_keep <- c()
  list_to_exclude <- c()
  for (i in 1:nrow(merged)){
    if ( (merged$V5[i] %in% c("A","C","G","T")==F) || (merged$V6[i] %in% c("A","C","G","T")==F) ){
      list_to_exclude <- c(list_to_exclude, paste(merged$V2[i]))
    }
    else{
      list_to_keep <- c(list_to_keep, paste(merged$V2[i]))
      positions_to_keep <- c(positions_to_keep, paste(merged$pos[i]))
    }
  }
  dup <- data.frame(table(positions_to_keep))
  positions_to_exclude <- dup$positions_to_keep[which(dup$Freq > 1)]
  list_to_exclude <- c(list_to_exclude, paste(list_to_keep[which(positions_to_keep %in% positions_to_exclude)]))
  write.table(list_to_exclude, file="list_to_exclude.txt", col.names = F, row.names = F, quote = F)
  rm(list = ls(all.names = TRUE))
}else{
  # just make an empty file
  write.table(loc_dup, file="list_to_exclude.txt", col.names = F, row.names = F, quote = F)
}

EOF
Rscript rem_loc_dup.R
# Step 4: exclude those SNPs
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar --exclude list_to_exclude.txt --make-bed --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude

# Step 4: If the next line outputs nothing, it means that we only have unique locations in the bim file. So the following should output nothing!!!
awk '{print $2}' hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude.bim | sort | uniq -d

# Step 5: update snp names to chr:pos if you want to
echo "# *-*-*-* update snp names to chr:pos"
awk '{print $2, $1 ":" $4}' hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude.bim > snps.txt
wc -l snps.txt
head snps.txt
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude --update-name snps.txt --make-bed --allow-extra-chr --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed
head hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed.bim

# Step 6: Make a set of the plink files that only has unrelateds (for training the models)
# Step 6: We remove 1st and 2nd relateds, while retaining as many people as possible
echo "# *-*-*-* Make a set of the plink files that only has unrelateds (for training the models)"
cat << 'EOF' >c_subset_unrelateds.R
library(bigsnpr)
library(dplyr)
# you need to give your path to the plink executable
plink2 <- "./PRS_software/plink2"
plink <- "./PRS_software/plink"
bedfile <- "hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed.bed"
# unlink("hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel.*")
# LD pruning is not recommended in KING
bedfile2 <- snp_plinkKINGQC(plink2, bedfile, thr.king = 2^-3.5, verbose = TRUE)
EOF
Rscript c_subset_unrelateds.R

# Step 6: this script will automatically output files of unrelateds with _norel appended to the file name:
# hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel.*

###### TRAINING OPTION 1 ######
# Make the training dataset. If you are using an external GWAS, you could just wrangle the sumstats into the format I describe below
###### TRAINING OPTION 2 ######
# But if you are running your own GWAS you can use those effect sizes (e.g. from plink --assoc or --linear etc) this code might help
# This is 100% necessary but I like to make a complete set of summary stats with allele frequency info
# Start by making a vcf from your plink files
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel --recode vcf --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel
# get the --freq info for the variants using vcftools
./vcftools --gzvcf hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel.vcf --freq --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel
# output file will be named hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel.frq
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel --assoc --allow-no-sex --out hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel
# output summary statistics will be named hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed_norel.assoc
# if you want a more traditional summary stats file 
# you can then wrangle and merge the *.frq and the *.assoc files by position
###### IMPORTANT TO CONSIDER ######
# Be careful to not train and validate (generate PRSs) within the same population
# That would lead to overfitting (if you ran associations and then computed genetic scores using the same set of individuals)
# A common workaround would be to run the assocaitions in 80% of the cohort
# And then vlaidate in 20% of the cohort
# And then perform crossvalidations for that 80% training 20% testing split
###### P VALUE THRESHOLDING ######
# If you want to subset your training dataset to a set of more significant SNPs do so now...
###### FORMATTING ######
# the resulting training effect sizes should be formatted like this...
# the first column matches the .bim SNP names
# the second column is the effect allele from the GWAS
# the third column is the log odds ratio, so the beta values 
# The file has no header. Let's name this file training_file.txt
# 1:152277622     T      0.5122
# 1:152280023     A      0.6773
# 1:152285076     CACTG  -1.0332
# 1:152285861     A      0.7526
###### MAKING THE PRSs ######
# Do some LDprunning, idk if this is the right amount of pruning.
# You should check the --indep parameters
# Note, for scoring we can use all individuals, not just the set of "unrelateds"
# But to measure PRS accuracy you would want to only compute AUC, R^2 etc using the unrelateds
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed --allow-no-sex --indep 50 5 5 --make-bed --out myplink_ldprune
plink --bfile hwe_geno_maf_rm_high_miss_SNPs_rm_hetFailures_remDupVar_exclude_renamed --extract myplink_ldprune.prune.in --make-bed --out plink_for_PRS
# Validate using plink --score and you will get an output file of PRSs
# The score is simply a weighted sum across SNPs 
# e.g. for a SNP take the number of an individual's specified effect alleles (0,1 or 2) and multiply by the weight (beta)
# then sum accross all SNPs. For a full example look at:
# https://zzz.bwh.harvard.edu/plink/profile.shtml
plink --bfile plink_for_PRS --score training_file.txt --out plink_score
# running the command above would generate a file plink_score.profile
# with one individual per row and the fields:
#      FID     Family ID
#      IID     Individual ID
#      PHENO   Phenotype for that
#      CNT     Number of non-missing SNPs used for scoring
#      CNT2    The number of named alleles
#      SCORE   Total score for that individual

