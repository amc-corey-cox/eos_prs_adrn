
library(ggplot2)
library(scales)
library(ggsignif)
library(pROC)
detach("package:pROC", unload=TRUE)
library(pROC)
library(fmsb)

#YOU ONLY NEED TO CHANGE THESE FIRST FEW THINGS and then just source it!

workingDirectory="~/Repositories/eos_prs_adrn/"
# "/wherever/your/files/are"

PRSfile="prs_file.txt"
# IID,true_phens,PRS
# 12341234,1,0.112676037532199
# 62324555,1,2.41159394913763
# 93845747,1,1.53975028442096

weightsFile="severity_file.txt"
# Normally I would have something like this:
# IID,Age,Sex,ancestry
# 12341234,64,Male,European
# 62324555,64,Female,European
# 93845747,42,Female,European

sev="EOS"
#the name of the column in weightsFile that you want to scale point colors with (also the name of the legend)
#set weightsFile=FALSE if you dont have a secondary file:
# weightsFile=FALSE

sevTreshold=F
#if you want to only keep individuals who have sev >= sevThreshold, else set =FALSE to keep everyone

outputPDF = "EA_ALL_EOS_PRS_vs_AD_Case_Control_violinPlot.pdf"
#whatever you want to name your pdf
# This script will generate a summary file titled: violinPlotter_summary.csv

plotTitle = "EA ALL -- adjusted EOS PRS vs AD -- Case/Control"

orLabel = TRUE
aucLabel = TRUE
pvalLabel = TRUE
plotLegend = TRUE
weightsLabel = TRUE
plotPointsUniform=FALSE
thePointColor="snow4"
# whatever you want to colour your points as. 
# This will be the NA color if points are missing, or the color of all pointsif plotPointsUniform==TRUE
plotPoints=TRUE
# if you want points to be in the plot
addToTitle="All AD"
#whatever you want the title to be
perturbFlag <- F
scaleFlag <- T
minMaxFlag <- F
upbdFlag <- F
PRSsoftware="LDpred"
# PRSsoftware="PLINK"

plotAUC = T
#set = TRUE if you want to generate an initial plot describe the AUC curve
weightsRegression = FALSE
#if you care to see the results of how well the weights you provided predict PRS

label_disty = .13
#The distance you want the AUC and OR brackets from min and max points on the plot
# A decent default is around 0.1

#               YOU SHOULDNT HAVE TO CHANGE ANYTHING BELOW THIS               #
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#               YOU SHOULDNT HAVE TO CHANGE ANYTHING BELOW THIS               #

setwd(workingDirectory)
# pdf(outputPDF, width = 8.5, height = 11 )

unlink("violinPlotter_summary.csv")
summaryTable <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(summaryTable) <- c("addToTitle","SevMeasure","OR Control vs Case", "OR LCI Control vs Case", "OR UCI Control vs Case", "AUC Control vs Case", "GLM pval Control vs Case", "sev Control Beta", "sev Case Beta","sev Control pval", "sev Case pval")
suppressWarnings(write.table(summaryTable, "violinPlotter_summary.csv", sep = ",", col.names = !file.exists("violinPlotter_summary.csv"), append = T, row.names = F))

# for (m in 1:length(PRSfileList)){
#   PRSfile <- PRSfileList[m]
#   addToTitle <- PRSfileList[m]

#load PRS file
if(PRSsoftware=="LDpred"){
  PRStable <- read.table(PRSfile,header=T,sep=",")
  colnames(PRStable) <- c("IID","PHENO","SCORE")
}else if(PRSsoftware=="XP-BLUP" || PRSsoftware=="PLINK"){
  PRStable <- read.table(PRSfile,header=T)
}else{
  print(paste0("Your PRSsoftware was: ",PRSsoftware))
  stop("ERROR: Your PRSsoftware must be one of LDpred, XP-BLUP, PLINK")
}

#load points
if(weightsFile!=FALSE){
  pointWeights <- read.table(weightsFile, header = T, stringsAsFactors = F, sep=",")
  pointWeights <- sapply(pointWeights, as.character)
  pointWeights[is.na(pointWeights)] <- ""
  pointWeights[pointWeights=="NA"] <- ""
  pointWeights[pointWeights==""] <- NA
  pointWeights <- as.data.frame(pointWeights)
}else if(weightsFile==FALSE){
  pointWeights <- PRStable
  pointWeights$SCORE <- NULL
  pointWeights$PHENO <- NULL
  pointWeights$weights <- c(1:length(pointWeights$weights))
  colnames(pointWeights) <- c("IID",sev)
  weightsRegression = FALSE
  plotLegend = FALSE
  weightsLabel = FALSE
  plotPointsUniform=TRUE
}

#merge the data together
LDpred_merge <- merge(PRStable, pointWeights, by=c("IID"))
LDpred_merge$Phen_num12 <- LDpred_merge$PHENO - 1
LDpred_merge$PHENO <- gsub('1', 'Control', LDpred_merge$PHENO)
LDpred_merge$PHENO <- gsub('2', 'Case', LDpred_merge$PHENO)
if(perturbFlag == T){
  perturb <- runif(n=nrow(LDpred_merge),min = -0.001,max = 0.001)
  LDpred_merge$SCORE <- LDpred_merge$SCORE + perturb
}

if(sevTreshold!=FALSE){
  LDpred_merge[[sev]] <- as.character(LDpred_merge[[sev]])
  LDpred_merge[[sev]] <- as.numeric(LDpred_merge[[sev]])
  LDpred_merge <- subset(LDpred_merge, LDpred_merge[[sev]]>=sevTreshold | is.na(LDpred_merge[[sev]])==TRUE)
}else{
  LDpred_merge[[sev]] <- as.character(LDpred_merge[[sev]])
  LDpred_merge[[sev]] <- as.numeric(LDpred_merge[[sev]])
}

#run binomial glm to find significance
if(scaleFlag == T){
  LDpred_merge$SCORE <- scale(LDpred_merge$SCORE)
}
# LDpred_merge$Phen_num12 <- as.numeric(LDpred_merge$PHENO)-1
# LDpred_merge$Phen_num12[LDpred_merge$Phen_num12==3] <- "2"
# LDpred_merge$Phen_num12 <- as.factor(LDpred_merge$Phen_num12)
# LDpred_merge$PHENO <- NA
# LDpred_merge$PHENO[LDpred_merge$Phen_num12==2] <- "Case"
# LDpred_merge$PHENO[LDpred_merge$Phen_num12==1] <- "Control"
LDpred_merge$PHENO <- as.factor(LDpred_merge$PHENO)
#case vs control
# LDpred_merge$SCORE <- -1*LDpred_merge$SCORE
# LDpred_merge$AgeSQ <- (LDpred_merge$Age)^2
caseControlGLM <- glm(Phen_num12 ~ SCORE, data = LDpred_merge, family = "binomial", na.action = na.omit)
# caseControlGLM <- glm(Phen_num12 ~ SCORE + Age + AgeSQ + Sex, data = LDpred_merge, family = "binomial", na.action = na.omit)
predpr <- predict(caseControlGLM,type=c("response"))
caseControlroccurve <- roc(LDpred_merge$Phen_num12[!is.na(LDpred_merge$Phen_num12)] ~ predpr)
caseControlroccurveCI <- roc(LDpred_merge$Phen_num12[!is.na(LDpred_merge$Phen_num12)] ~ predpr, ci=T)

auc(case ~ age, data=infert)   #Compute AUC for predicting case with the variable age
mod1<-glm(case ~ age + parity, data=infert, family="binomial")  #Logistic regression model
auc(case ~ predict(mod1), data=infert)  #Compute AUC for predicting case with your model

if (plotAUC==T){
  caseControlplot <- plot(caseControlroccurve, main=paste("Case vs Control AUC =", round(caseControlroccurve$auc, 3)))
}
caseControlp <- formatC(coef(summary(caseControlGLM))[,4][2], format = "e", digits = 0)
caseControlOR <- exp(cbind("Odds ratio" = coef(caseControlGLM), confint.default(caseControlGLM, level = 0.95)))
if (as.numeric(strsplit(caseControlp,"")[[1]][1]) == 1){
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
}else{
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p < ", 10*as.numeric(caseControlp)/as.numeric(strsplit(caseControlp,"")[[1]][1]))
}  
if (as.numeric(caseControlp) >= 0.001){
  caseControlpthresh <- paste0("OR [95% CI] = ", format(round(as.numeric(caseControlOR[2,1]), 2), nsmall = 2), " [", format(round(as.numeric(caseControlOR[2,2]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlOR[2,3]), 2),nsmal = 2),"], ","p = ", formatC(as.numeric(caseControlp), format = "g"))
}
caseControlORsummary <- paste0("AUC [95% CI] = ", format(round(as.numeric(caseControlroccurve$auc), 2), nsmall = 2), " [", format(round(as.numeric(caseControlroccurveCI$ci[1]), 2),nsmal = 2), ", ", format(round(as.numeric(caseControlroccurveCI$ci[3]), 2),nsmal = 2),"]")

if (weightsRegression == FALSE || sum(is.na(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"]))==length(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"])){
  ControlsevlmB <- "NA"
  Controlsevlmp <- "NA"
}else{
  Controlsevlm <- lm(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"] ~ LDpred_merge$SCORE[LDpred_merge$PHENO == "Control"])
  ControlsevlmB <- format(summary(Controlsevlm)$coefficients[2,1],digits=2)
  Controlsevlmp <- format(summary(Controlsevlm)$coefficients[2,4],digits=2)
}
if (weightsRegression == FALSE || sum(is.na(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"]))==length(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"])){
  CasesevlmB <- "NA"
  Casesevlmp <- "NA"
}else{
  Casesevlm <- lm(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"] ~ LDpred_merge$SCORE[LDpred_merge$PHENO == "Case"])
  CasesevlmB <- format(summary(Casesevlm)$coefficients[2,1],digits=2)
  Casesevlmp <- format(summary(Casesevlm)$coefficients[2,4],digits=2)
}

if (weightsRegression == TRUE){
  Bothlm <- lm(LDpred_merge[[sev]] ~ LDpred_merge$SCORE)
  BothlmB <- format(summary(Casesevlm)$coefficients[2,1],digits=2)
  Bothlmp <- format(summary(Casesevlm)$coefficients[2,4],digits=2)
  sink("weights_LM_Summary.txt")
  print("This is a summary of the LM for how well the weights you provided predict the stanardized PRS")
  print(summary(Bothlm))
  sink()  # returns output to the console
}

#set up the data, and the order in which we want to do the violin plot
LDpred_merge$PHENO <- factor(LDpred_merge$PHENO, levels = c("Control","Case"),ordered = T)

minSev <- min(LDpred_merge[[sev]],na.rm=T)
maxSev <- max(LDpred_merge[[sev]],na.rm=T)

if(minMaxFlag == T){
  a <- -3.380482
  b <- 5
  minPRS <- a-(abs(b-a))*(label_disty+.002)
  maxPRS <- b + (abs(b-a))*(label_disty+.008)
  botLabLoc <- a-(abs(b-a))*(label_disty)
  topLabLoc <- b + (abs(b-a))*(label_disty*.82)
}else{
  minPRS <- min(LDpred_merge$SCORE)-(abs(max(LDpred_merge$SCORE)-min(LDpred_merge$SCORE)))*(label_disty+.002)
  maxPRS <- max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE)-min(LDpred_merge$SCORE)))*(label_disty+.008)
  botLabLoc <- min(LDpred_merge$SCORE)-(abs(max(LDpred_merge$SCORE)-min(LDpred_merge$SCORE)))*(label_disty)
  topLabLoc <- max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE)-min(LDpred_merge$SCORE)))*(label_disty*.82)
}

LDpredPlot <-
  ggplot(data = LDpred_merge,aes(x = PHENO, y = SCORE))+
  scale_fill_viridis_d(option = "D")+
  theme_dark(base_size = 8)+
  geom_violin(fill = "gray70",alpha=0.4, position = position_dodge(width = .5),size=1,color="gray22",width=.5,lwd=.2) +
  geom_boxplot(fill = "gray95",notch = F, shape=21, outlier.size = -1, color="gray32",lwd=.5, alpha = .75)+
  theme(plot.title = element_text(size=8))+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_text(size=8))+
  theme(axis.text.x = element_text(colour = "black",size=8))+
  theme(axis.text.y = element_text(colour = "black",size=7.6))+
  theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"))+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color='gray80'))+
  # labs(title = addToTitle)+
  ylab("Standardized PRS")+
  scale_y_continuous(breaks = seq(-100,100, by=2), limits = c(minPRS,maxPRS)) +
  labs(title = plotTitle)

if(orLabel == TRUE && aucLabel == TRUE && pvalLabel == TRUE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlORsummary, color="black", y_position = topLabLoc,tip_length=.03)+
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03)
}else if(orLabel == TRUE && aucLabel == FALSE && pvalLabel == TRUE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlORsummary, color="black", y_position = topLabLoc,tip_length=.03)+
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=strsplit(caseControlpthresh,", ")[[1]][1], color="black", y_position = botLabLoc,tip_length=-.03)
}else if(orLabel == TRUE && aucLabel == TRUE && pvalLabel == FALSE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlORsummary, color="black", y_position = topLabLoc,tip_length=.03)+
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=strsplit(caseControlpthresh,", ")[[1]][2], color="black", y_position = botLabLoc,tip_length=-.03)
}else if(orLabel == TRUE && aucLabel == FALSE && pvalLabel == FALSE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlORsummary, color="black", y_position = topLabLoc,tip_length=.03)
}else if(orLabel == FALSE && aucLabel == TRUE && pvalLabel == TRUE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=caseControlpthresh, color="black", y_position = botLabLoc,tip_length=-.03)
}else if(orLabel == FALSE && aucLabel == TRUE && pvalLabel == FALSE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=strsplit(caseControlpthresh,", ")[[1]][2], color="black", y_position = botLabLoc,tip_length=-.03)
}else if(orLabel == FALSE && aucLabel == FALSE && pvalLabel == TRUE){
  LDpredPlot <- LDpredPlot +
    geom_signif(textsize = 2.25, comparisons = list(c("Case", "Control")), annotations=strsplit(caseControlpthresh,", ")[[1]][1], color="black", y_position = botLabLoc,tip_length=-.03)
}

if(plotPointsUniform == TRUE && plotPoints == TRUE){
  plotLegend = FALSE
  LDpredPlot <- LDpredPlot +
    # geom_point(aes(col = LDpred_merge[[sev]]), size=.05, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = F)+
    # geom_point(aes(col = LDpred_merge[[sev]]), size=1, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = F)+
    scale_colour_gradient(low = thePointColor, high = thePointColor, na.value=thePointColor)+
    labs(colour=sev)
}else{
  if(plotLegend == FALSE && plotPoints == TRUE){
    LDpredPlot <- LDpredPlot +
      geom_point( aes(col = LDpred_merge[[sev]]), size=.05, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = F)+
      # geom_point( aes(col = LDpred_merge[[sev]]), size=1, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = F)+
      scale_colour_gradientn(name = paste0("", sev),colours = rev(c("#D7191C","#D7191C","#FDAE61","#ABDDA4","#2B83BA","#2B83BA")),values = rescale(c(-0.03780071,-0.01890036,0.422549,1.112655,1.624349,1.868488), to = c(0, 1)), na.value=thePointColor,limits=c(minSev,maxSev))+
      labs(colour=sev)
  }else if(plotLegend == TRUE && plotPoints == TRUE){
    LDpredPlot <- LDpredPlot +
      geom_point( aes(col = LDpred_merge[[sev]]), size=.05, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = T)+
      # geom_point( aes(col = LDpred_merge[[sev]]), size=1, position = position_jitterdodge(seed = 1,jitter.width = .5), alpha=.3, show.legend = T)+
      scale_colour_gradientn(name = paste0("", sev),colours = rev(c("#D7191C","#D7191C","#FDAE61","#ABDDA4","#2B83BA","#2B83BA")),values = rescale(c(-0.03780071,-0.01890036,0.422549,1.112655,1.624349,1.868488), to = c(0, 1)), na.value=thePointColor,limits=c(minSev,maxSev))+
      labs(colour=sev)+
      theme(legend.title = element_text(size = 7), 
            legend.text = element_text(size = 7))+
      theme(legend.key.size = unit(0.25, "cm"))
  }
}

if(weightsLabel == TRUE && weightsRegression == TRUE){
  LDpredPlot <- LDpredPlot +
    xlab(bquote(paste(.("Control "),beta,.(" = "),.(ControlsevlmB),.(" p = "),.(Controlsevlmp),.("  |  Case "),beta,.(" = "),.(CasesevlmB),.(" p = "),.(Casesevlmp), sep = "")))
}else{
  LDpredPlot <- LDpredPlot +
    xlab("")
}

V <- ggplot_build(LDpredPlot)

pdf(outputPDF,width=4,height=3.5)
print(LDpredPlot)
# print(caseControlplot)
dev.off()
# save(LDpred_merge,LDpredPlot,sev, file="V.Rdata",PRSfile)
# dev.off()

summaryTable <- data.frame(matrix(ncol = 11, nrow = 1))
colnames(summaryTable) <- c("addToTitle","SevMeasure","OR Control vs Case", "OR LCI Control vs Case", "OR UCI Control vs Case", "AUC Control vs Case", "GLM pval Control vs Case", "sev Control Beta", "sev Case Beta","sev Control pval", "sev Case pval")
summaryTable[1,] <- c(addToTitle,
                      sev,
                      format(round(as.numeric(caseControlOR[2,1]), 2)), 
                      format(round(as.numeric(caseControlOR[2,2]), 2)), 
                      format(round(as.numeric(caseControlOR[2,3]), 2)), 
                      round(caseControlroccurve$auc, 2),
                      as.numeric(caseControlp),
                      ControlsevlmB,
                      CasesevlmB,
                      Controlsevlmp, 
                      Casesevlmp)

write.table(summaryTable, "violinPlotter_summary.csv", sep = ",", col.names = !file.exists("violinPlotter_summary.csv"), append = T, row.names = F)
# }

# Get the Nagelkerke and McFadden R2 estimates
mod <- glm(Phen_num12 ~ SCORE, data = LDpred_merge, family = "binomial", na.action = na.omit)
nullmod <- glm(Phen_num12~1, data = LDpred_merge, family="binomial", na.action = na.omit)
1-logLik(mod)/logLik(nullmod)
1-mod$deviance/mod$null.deviance
NagelkerkeR2(mod)

print(LDpredPlot)

