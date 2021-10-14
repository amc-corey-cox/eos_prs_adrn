gg_prs_prep_data <- function(prs_data, score, pheno, ...) {
  LDpred_merge <- prs_data %>%
    mutate(
      PHENO = !!ensym(pheno),
      SCORE = !!ensym(score),
      PHENO_NUM_01 = PHENO - 1,
      PHENO = recode_factor(PHENO, `1` = "Control", `2` = "Case", .ordered = TRUE)
    ) %>%
    select(PHENO, SCORE, PHENO_NUM_01, !!! ensyms(...))
}

gg_prs_auc <- function(prs_data) {
  caseControlGLM <- glm(PHENO_NUM_01 ~ SCORE, data = prs_data, family = "binomial", na.action = na.omit)
  
  predpr <- predict(caseControlGLM, type = c("response"))
  prs_clean <- prs_data %>% filter(! is.na(PHENO_NUM_01))
  caseControlroccurve <- roc(prs_clean$PHENO_NUM_01 ~ predpr)
  caseControlroccurveCI <- roc(prs_clean$PHENO_NUM_01 ~ predpr, ci = TRUE)
  
  plot(caseControlroccurve, main = paste("Case vs Control AUC =", round(caseControlroccurve$auc, 3)))
}

gg_prs_violin <- function(prs_data, sev, plot_title, threshold = -Inf, use_scale = TRUE, perturb = FALSE, label_dist_y = 0.13) {
  # minMaxFlag <- FALSE
  
  perturb_vals  <- runif(n = nrow(prs_data), min = -0.001, max = 0.001)
  
  LDpred_merge <- prs_data %>%
    mutate(SCORE = if(use_scale) scale(SCORE) else SCORE,
           SCORE = if(perturb) SCORE + perturb_vals else SCORE) %>%
    filter(!! ensym(sev) > threshold)
  
  # if (minMaxFlag == TRUE) {
  #   a <- -3.380482
  #   b <- 5
  #   minPRS <- a - (abs(b - a)) * (label_dist_y + .002)
  #   maxPRS <- b + (abs(b - a)) * (label_dist_y + .008)
  #   botLabLoc <- a - (abs(b - a)) * (label_dist_y)
  #   topLabLoc <- b + (abs(b - a)) * (label_dist_y * .82)
  # } else{
  #   minPRS <-
  #     min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y + .002)
  #   maxPRS <-
  #     max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y + .008)
  #   botLabLoc <-
  #     min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y)
  #   topLabLoc <-
  #     max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y * .82)
  # }
  
  minPRS <-
    min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
    (label_dist_y + .002)
  maxPRS <-
    max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
    (label_dist_y + .008)
  
  LDpredPlot <-
    ggplot(data = LDpred_merge, aes(x = PHENO, y = SCORE)) +
    scale_fill_viridis_d(option = "D") +
    theme_dark(base_size = 8) +
    geom_violin(fill = "gray70", alpha = 0.4, position = position_dodge(width = .5), size = 1, color = "gray22", width = .5, lwd = .2) +
    geom_boxplot(fill = "gray95", notch = F, shape = 21, outlier.size = -1, color = "gray32", lwd = .5, alpha = .75) +
    theme(plot.title = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.x = element_text(size = 8)) +
    theme(axis.text.x = element_text(colour = "black", size = 8)) +
    theme(axis.text.y = element_text(colour = "black", size = 7.6)) +
    theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
    theme(panel.grid.major.x = element_blank()) +
    theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'gray80')) +
    ylab("Standardized PRS") +
    scale_y_continuous(breaks = seq(-100, 100, by = 2), limits = c(minPRS, maxPRS))
  
  LDpredPlot
}

gg_prs_violin_labels <- function (LDpredPlot, LDpred_merge, sev, plot_dist_y = 0.13) {
  
  # Flags and configuration 
  sevThreshold = FALSE
  #if you want to only keep individuals who have sev >= sevThreshold, else set =FALSE to keep everyone
  
  orLabel = TRUE
  aucLabel = TRUE
  pvalLabel = TRUE
  plotLegend = TRUE
  weightsLabel = TRUE
  plotPointsUniform = FALSE
  thePointColor = "snow4"
  # whatever you want to colour your points as.
  # This will be the NA color if points are missing, or the color of all pointsif plotPointsUniform==TRUE
  plotPoints = TRUE
  # if you want points to be in the plot
  perturbFlag <- FALSE
  scaleFlag <- TRUE
  minMaxFlag <- FALSE
  upbdFlag <- FALSE
  
  # plotAUC = TRUE
  #set = TRUE if you want to generate an initial plot describe the AUC curve
  weightsRegression = FALSE
  #if you care to see the results of how well the weights you provided predict PRS
  
  label_dist_y = .13
  #The distance you want the AUC and OR brackets from min and max points on the plot
  # A decent default is around 0.1
  
  ### Code begins here
  #run binomial glm to find significance
  
  botLabLoc <-
    min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
    (label_dist_y)
  topLabLoc <-
    max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
    (label_dist_y * .82)
  
  caseControlGLM <- glm(PHENO_NUM_01 ~ SCORE, data = LDpred_merge, family = "binomial", na.action = na.omit)
  
  predpr <- predict(caseControlGLM, type = c("response"))
  LDpred_clean <- LDpred_merge %>% filter(! is.na(PHENO_NUM_01))
  caseControlroccurve <- roc(LDpred_clean$PHENO_NUM_01 ~ predpr)
  caseControlroccurveCI <- roc(LDpred_clean$PHENO_NUM_01 ~ predpr, ci = TRUE)
  
  caseControlp <-
    formatC(coef(summary(caseControlGLM))[, 4][2], format = "e", digits = 0)
  caseControlOR <-
    exp(cbind(
      "Odds ratio" = coef(caseControlGLM),
      confint.default(caseControlGLM, level = 0.95)
    ))
  if (as.numeric(strsplit(caseControlp, "")[[1]][1]) == 1) {
    caseControlpthresh <-
      paste0(
        "OR [95% CI] = ",
        format(round(as.numeric(caseControlOR[2, 1]), 2), nsmall = 2),
        " [",
        format(round(as.numeric(caseControlOR[2, 2]), 2), nsmal = 2),
        ", ",
        format(round(as.numeric(caseControlOR[2, 3]), 2), nsmal = 2),
        "], ",
        "p < ",
        as.numeric(caseControlp) / as.numeric(strsplit(caseControlp, "")[[1]][1])
      )
  } else{
    caseControlpthresh <-
      paste0(
        "OR [95% CI] = ",
        format(round(as.numeric(caseControlOR[2, 1]), 2), nsmall = 2),
        " [",
        format(round(as.numeric(caseControlOR[2, 2]), 2), nsmal = 2),
        ", ",
        format(round(as.numeric(caseControlOR[2, 3]), 2), nsmal = 2),
        "], ",
        "p < ",
        10 * as.numeric(caseControlp) / as.numeric(strsplit(caseControlp, "")[[1]][1])
      )
  }
  if (as.numeric(caseControlp) >= 0.001) {
    caseControlpthresh <-
      paste0(
        "OR [95% CI] = ",
        format(round(as.numeric(caseControlOR[2, 1]), 2), nsmall = 2),
        " [",
        format(round(as.numeric(caseControlOR[2, 2]), 2), nsmal = 2),
        ", ",
        format(round(as.numeric(caseControlOR[2, 3]), 2), nsmal = 2),
        "], ",
        "p = ",
        formatC(as.numeric(caseControlp), format = "g")
      )
  }
  caseControlORsummary <-
    paste0(
      "AUC [95% CI] = ",
      format(round(
        as.numeric(caseControlroccurve$auc), 2
      ), nsmall = 2),
      " [",
      format(round(
        as.numeric(caseControlroccurveCI$ci[1]), 2
      ), nsmal = 2),
      ", ",
      format(round(
        as.numeric(caseControlroccurveCI$ci[3]), 2
      ), nsmal = 2),
      "]"
    )
  
  if (weightsRegression == FALSE ||
      sum(is.na(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"])) == length(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"])) {
    ControlsevlmB <- "NA"
    Controlsevlmp <- "NA"
  } else{
    Controlsevlm <-
      lm(LDpred_merge[[sev]][LDpred_merge$PHENO == "Control"] ~ LDpred_merge$SCORE[LDpred_merge$PHENO == "Control"])
    ControlsevlmB <-
      format(summary(Controlsevlm)$coefficients[2, 1], digits = 2)
    Controlsevlmp <-
      format(summary(Controlsevlm)$coefficients[2, 4], digits = 2)
  }
  if (weightsRegression == FALSE ||
      sum(is.na(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"])) == length(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"])) {
    CasesevlmB <- "NA"
    Casesevlmp <- "NA"
  } else{
    Casesevlm <-
      lm(LDpred_merge[[sev]][LDpred_merge$PHENO == "Case"] ~ LDpred_merge$SCORE[LDpred_merge$PHENO == "Case"])
    CasesevlmB <-
      format(summary(Casesevlm)$coefficients[2, 1], digits = 2)
    Casesevlmp <-
      format(summary(Casesevlm)$coefficients[2, 4], digits = 2)
  }
  
  if (weightsRegression == TRUE) {
    Bothlm <- lm(LDpred_merge[[sev]] ~ LDpred_merge$SCORE)
    BothlmB <- format(summary(Casesevlm)$coefficients[2, 1], digits = 2)
    Bothlmp <- format(summary(Casesevlm)$coefficients[2, 4], digits = 2)
    sink("weights_LM_Summary.txt")
    print(
      "This is a summary of the LM for how well the weights you provided predict the stanardized PRS"
    )
    print(summary(Bothlm))
    sink()  # returns output to the console
  }
  
  #set up the data, and the order in which we want to do the violin plot
  # LDpred_merge$PHENO <-
  #   factor(LDpred_merge$PHENO,
  #          levels = c("Control", "Case"),
  #          ordered = T)
  
  minSev <- min(LDpred_merge[[sev]], na.rm = T)
  maxSev <- max(LDpred_merge[[sev]], na.rm = T)
  
  # if (minMaxFlag == TRUE) {
  #   a <- -3.380482
  #   b <- 5
  #   minPRS <- a - (abs(b - a)) * (label_dist_y + .002)
  #   maxPRS <- b + (abs(b - a)) * (label_dist_y + .008)
  #   botLabLoc <- a - (abs(b - a)) * (label_dist_y)
  #   topLabLoc <- b + (abs(b - a)) * (label_dist_y * .82)
  # } else{
  #   minPRS <-
  #     min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y + .002)
  #   maxPRS <-
  #     max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y + .008)
  #   botLabLoc <-
  #     min(LDpred_merge$SCORE) - (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y)
  #   topLabLoc <-
  #     max(LDpred_merge$SCORE) + (abs(max(LDpred_merge$SCORE) - min(LDpred_merge$SCORE))) *
  #     (label_dist_y * .82)
  # }
  # 
  # LDpredPlot <-
  #   ggplot(data = LDpred_merge, aes(x = PHENO, y = SCORE)) +
  #   scale_fill_viridis_d(option = "D") +
  #   theme_dark(base_size = 8) +
  #   geom_violin(fill = "gray70", alpha = 0.4, position = position_dodge(width = .5), size = 1, color = "gray22", width = .5, lwd = .2) +
  #   geom_boxplot(fill = "gray95", notch = F, shape = 21, outlier.size = -1, color = "gray32", lwd = .5, alpha = .75) +
  #   theme(plot.title = element_text(size = 8)) + theme(plot.title = element_text(hjust = 0.5)) +
  #   theme(axis.title.x = element_text(size = 8)) +
  #   theme(axis.text.x = element_text(colour = "black", size = 8)) +
  #   theme(axis.text.y = element_text(colour = "black", size = 7.6)) +
  #   theme(axis.line = element_line(colour = "black", size = 0.5, linetype = "solid")) +
  #   theme(panel.grid.major.x = element_blank()) +
  #   theme(panel.background = element_rect(fill = 'white'), panel.grid = element_line(color = 'gray80')) +
  #   ylab("Standardized PRS") +
  #   scale_y_continuous(breaks = seq(-100, 100, by = 2), limits = c(minPRS, maxPRS))
  
  # print(LDpredPlot)
  
  if (orLabel == TRUE && aucLabel == TRUE && pvalLabel == TRUE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlORsummary,
        color = "black",
        y_position = topLabLoc,
        tip_length = .03
      ) +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlpthresh,
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  } else if (orLabel == TRUE &&
             aucLabel == FALSE && pvalLabel == TRUE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlORsummary,
        color = "black",
        y_position = topLabLoc,
        tip_length = .03
      ) +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = strsplit(caseControlpthresh, ", ")[[1]][1],
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  } else if (orLabel == TRUE &&
             aucLabel == TRUE && pvalLabel == FALSE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlORsummary,
        color = "black",
        y_position = topLabLoc,
        tip_length = .03
      ) +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = strsplit(caseControlpthresh, ", ")[[1]][2],
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  } else if (orLabel == TRUE &&
             aucLabel == FALSE && pvalLabel == FALSE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlORsummary,
        color = "black",
        y_position = topLabLoc,
        tip_length = .03
      )
  } else if (orLabel == FALSE &&
             aucLabel == TRUE && pvalLabel == TRUE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = caseControlpthresh,
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  } else if (orLabel == FALSE &&
             aucLabel == TRUE && pvalLabel == FALSE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = strsplit(caseControlpthresh, ", ")[[1]][2],
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  } else if (orLabel == FALSE &&
             aucLabel == FALSE && pvalLabel == TRUE) {
    LDpredPlot <- LDpredPlot +
      geom_signif(
        textsize = 2.25,
        comparisons = list(c("Case", "Control")),
        annotations = strsplit(caseControlpthresh, ", ")[[1]][1],
        color = "black",
        y_position = botLabLoc,
        tip_length = -.03
      )
  }
  
  if (plotPointsUniform == TRUE && plotPoints == TRUE) {
    plotLegend = FALSE
    LDpredPlot <- LDpredPlot +
      geom_point(
        aes(col = LDpred_merge[[sev]]),
        size = .05,
        position = position_jitterdodge(seed = 1, jitter.width = .5),
        alpha = .3,
        show.legend = F
      ) +
      scale_colour_gradient(low = thePointColor,
                            high = thePointColor,
                            na.value = thePointColor) +
      labs(colour = sev)
  } else{
    if (plotLegend == FALSE && plotPoints == TRUE) {
      LDpredPlot <- LDpredPlot +
        geom_point(
          aes(col = LDpred_merge[[sev]]),
          size = .05,
          position = position_jitterdodge(seed = 1, jitter.width = .5),
          alpha = .3,
          show.legend = F
        ) +
        scale_colour_gradientn(
          name = paste0("", sev),
          colours = rev(
            c(
              "#D7191C",
              "#D7191C",
              "#FDAE61",
              "#ABDDA4",
              "#2B83BA",
              "#2B83BA"
            )
          ),
          values = rescale(
            c(
              -0.03780071,
              -0.01890036,
              0.422549,
              1.112655,
              1.624349,
              1.868488
            ),
            to = c(0, 1)
          ),
          na.value = thePointColor,
          limits = c(minSev, maxSev)
        ) +
        labs(colour = sev)
    } else if (plotLegend == TRUE && plotPoints == TRUE) {
      LDpredPlot <- LDpredPlot +
        geom_point(
          aes(col = LDpred_merge[[sev]]),
          size = .05,
          position = position_jitterdodge(seed = 1, jitter.width = .5),
          alpha = .3,
          show.legend = T
        ) +
        scale_colour_gradientn(
          name = paste0("", sev),
          colours = rev(
            c(
              "#D7191C",
              "#D7191C",
              "#FDAE61",
              "#ABDDA4",
              "#2B83BA",
              "#2B83BA"
            )
          ),
          values = rescale(
            c(
              -0.03780071,
              -0.01890036,
              0.422549,
              1.112655,
              1.624349,
              1.868488
            ),
            to = c(0, 1)
          ),
          na.value = thePointColor,
          limits = c(minSev, maxSev)
        ) +
        labs(colour = sev) +
        theme(legend.title = element_text(size = 7),
              legend.text = element_text(size = 7)) +
        theme(legend.key.size = unit(0.25, "cm"))
    }
  }

  if (weightsLabel == TRUE && weightsRegression == TRUE) {
    LDpredPlot <- LDpredPlot +
      xlab(bquote(paste(.("Control "), beta, .(" = "), .(ControlsevlmB), .(" p = "), .(Controlsevlmp), .("  |  Case "), beta, .(" = "), .(CasesevlmB), .(" p = "), .(Casesevlmp), sep = "")))
  } else{
    LDpredPlot <- LDpredPlot + xlab("")
  }
  
  # V <- ggplot_build(LDpredPlot)
  
  # summaryTable <- data.frame(matrix(ncol = 11, nrow = 1))
  # colnames(summaryTable) <-
  #   c(
  #     "addToTitle",
  #     "SevMeasure",
  #     "OR Control vs Case",
  #     "OR LCI Control vs Case",
  #     "OR UCI Control vs Case",
  #     "AUC Control vs Case",
  #     "GLM pval Control vs Case",
  #     "sev Control Beta",
  #     "sev Case Beta",
  #     "sev Control pval",
  #     "sev Case pval"
  #   )
  # summaryTable[1, ] <- c(
  #   addToTitle,
  #   sev,
  #   format(round(as.numeric(caseControlOR[2, 1]), 2)),
  #   format(round(as.numeric(caseControlOR[2, 2]), 2)),
  #   format(round(as.numeric(caseControlOR[2, 3]), 2)),
  #   round(caseControlroccurve$auc, 2),
  #   as.numeric(caseControlp),
  #   ControlsevlmB,
  #   CasesevlmB,
  #   Controlsevlmp,
  #   Casesevlmp
  # )
  
  # Get the Nagelkerke and McFadden R2 estimates
  mod <-
    glm(
      PHENO_NUM_01 ~ SCORE,
      data = LDpred_merge,
      family = "binomial",
      na.action = na.omit
    )
  nullmod <-
    glm(
      PHENO_NUM_01 ~ 1,
      data = LDpred_merge,
      family = "binomial",
      na.action = na.omit
    )
  1 - logLik(mod) / logLik(nullmod)
  1 - mod$deviance / mod$null.deviance
  NagelkerkeR2(mod)
  
  LDpredPlot
}
