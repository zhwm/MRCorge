#' Perform Mendelian Randomization based on the core gene hypothesis for polygenic exposures
#'
#' This function performs Mendelian Randomization (MR) analyses based on the core gene hypothesis for polygenic exposures.
#' MR Corge identifies a small number of putative core instruments that are more likely to affect genes with a direct biological role 
#' in an exposure and obtains causal effect estimates based on these instruments, thereby reducing the risk of horizontal pleiotropy.
#'
#' @param harmonized_data Data frame containing the harmonized data from the TwoSampleMR package.
#' @param rank Character, method to rank instruments, one of 'beta' (absolute per-allele effect size), 'h2' (per-variant heritability), 'pval' (significance), or 'zscore' (sample size-normalized z-score). Default is 'beta'.
#' @param K Integer, number of instrument groups. Default is 5.
#' @param method_list Character vector, list of MR methods to apply. Default is c("mr_ivw", "mr_weighted_mode", "mr_weighted_median").
#' @param seed Integer, seed for random number generation. Default is 1.
#' @return A list containing:
#' \describe{
#'   \item{mrcorge}{The main MR results.}
#'   \item{result}{All MR results including group-specific and cumulative estimates.}
#'   \item{instrument}{Data frame of instruments with their respective groups.}
#'   \item{Q}{Q statistics for heterogeneity.}
#'   \item{I2}{I2 statistics for heterogeneity.}
#'   \item{het_P}{P-values for heterogeneity tests.}
#' }
#' @examples
#' # Example usage:
#' library(MRCorge)
#' data("HDL_CAD")
#' res <- mrcorge(HDL_CAD, rank='beta', K=20, method_list = c("mr_ivw"))
#' @export
mrcorge <- function(harmonized_data, rank = 'beta', K = 5, method_list = c("mr_ivw","mr_weighted_mode","mr_weighted_median"), seed = 1) {
  set.seed(seed)
  library(ggplot2)
  library(TwoSampleMR)
  harmonized_data <- subset(harmonized_data, mr_keep & !is.na(se.exposure) & !is.na(eaf.exposure))
  if (nrow(harmonized_data)<=K) {
    stop("K should be smaller than the number of candidate instruments in harmonized_data!")
  }
  # Rank instruments based on the chosen metric
  rank_values <- switch(rank,
                        'beta' = abs(harmonized_data$beta.exposure),
                        'h2' = (harmonized_data$beta.exposure / harmonized_data$se.exposure)^2 / 
                          ((harmonized_data$beta.exposure / harmonized_data$se.exposure)^2 + harmonized_data$samplesize.exposure - 2),
                        'pval' = abs(harmonized_data$beta.exposure / harmonized_data$se.exposure),
                        'zscore' = abs(harmonized_data$beta.exposure / harmonized_data$se.exposure / sqrt(harmonized_data$samplesize.exposure)),
                        stop("Please choose one of the ranking methods from c('beta', 'h2', 'pval', 'zscore')!"))
  
  harmonized_data$instrument_group <- K + 1 - as.numeric(cut_number(rank_values + rnorm(length(rank_values), sd = 1e-10), K))
  
  instrument_group <- with(harmonized_data, data.frame(
    instrument = SNP,
    beta.exposure = beta.exposure,
    eaf.exposure = eaf.exposure,
    rank_value = rank_values,
    instrument_group = instrument_group,
    chr = if("chr.exposure" %in% colnames(harmonized_data)) chr.exposure else NA,
    pos = if("pos.exposure" %in% colnames(harmonized_data)) pos.exposure else NA
  ))
  
  all_result <- data.frame(id.exposure = character(), id.outcome = character(), outcome = character(), exposure = character(), 
                           method = character(), nsnp = numeric(), b = numeric(), se = numeric(), pval = numeric(), 
                           Fstat = numeric(), group = numeric(), category = character())
  
  process_group <- function(data, group, category) {
    result_group <- mr(data, method_list = method_list)
    if ("mr_egger_regression" %in% method_list) {
      result_group_egger_intercept <- mr_egger_regression(data$beta.exposure, data$beta.outcome, data$se.exposure, data$se.outcome)
      result_group <- rbind(result_group, data.frame(
        id.exposure = result_group$id.exposure[1],
        id.outcome = result_group$id.outcome[1],
        outcome = result_group$outcome[1],
        exposure = result_group$exposure[1],
        method = "MR Egger intercept",
        nsnp = result_group_egger_intercept$nsnp,
        b = result_group_egger_intercept$b_i,
        se = result_group_egger_intercept$se_i,
        pval = result_group_egger_intercept$pval_i
      ))
    }
    R2 <- sum((data$beta.exposure / data$se.exposure)^2 / ((data$beta.exposure / data$se.exposure)^2 + data$samplesize.exposure - 2))
    result_group$Fstat <- (mean(data$samplesize.exposure) - nrow(data) - 1) / nrow(data) * R2 / (1 - R2)
    result_group$group <- group
    result_group$category <- category
    result_group
  }
  
  for (group in 1:K) {
    group_data <- harmonized_data[harmonized_data$instrument_group == group, ]
    group_res <- process_group(group_data, group, "Group-specific")
    if (group == 1) {
      cat("MR CORGE estimates: \n")
      print(group_res)
    } else {
      cat("obtaining group-specific estimates for instrument group", group, "/", K, "...", collapse = "\n")
    }
    all_result <- rbind(all_result, group_res)
  }
  
  for (group in 1:K) {
    cat("obtaining cumulative estimates for instrument group", group, "/", K, "...", collapse = "\n")
    group_data <- harmonized_data[harmonized_data$instrument_group <= group, ]
    all_result <- rbind(all_result, process_group(group_data, group, "Cumulative"))
  }
  
  Q <- sapply(unique(all_result$method[all_result$method != "MR Egger intercept"]), function(m) 
    Qstatistic(all_result[all_result$method == m & all_result$category == "Group-specific",]$b,
               all_result[all_result$method == m & all_result$category == "Group-specific",]$se,
               all_result[all_result$method == m & all_result$category == "Cumulative" & all_result$group == K,]$b))
  
  I2 <- pmax((Q - (K - 1)) / Q, 0)
  het_P <- pchisq(Q, df = K - 1, lower.tail = FALSE)
  
  list(mrcorge = all_result[all_result$category == "Group-specific" & all_result$group == 1, ],
       result = all_result,
       instrument = instrument_group,
       Q = Q,
       I2 = I2,
       het_P = het_P)
}

#' Calculate Q Statistic for Heterogeneity
#'
#' This function calculates the Q statistic for heterogeneity in Mendelian Randomization analysis.
#'
#' @param b_group Numeric vector, effect sizes for each group.
#' @param se_group Numeric vector, standard errors for each group.
#' @param b_overall Numeric, overall effect size.
#' @return Numeric, Q statistic for heterogeneity.
#' @export
Qstatistic <- function(b_group, se_group, b_overall) {
  sum(1 / se_group^2 * (b_group - b_overall)^2)
}

#' Plot MR CORGE Results
#'
#' This function plots the results of MR CORGE.
#'
#' @param mrcorge_result List, result from the mrcorge function.
#' @param scale Character, scale for the plot, either "linear" or "exp". Default is "linear".
#' @return ggplot object, the plot of MR sensitivity analysis results.
#' @examples
#' # Example usage:
#' library(MRCorge)
#' data("HDL_CAD")
#' res <- mrcorge(HDL_CAD, rank='beta', K=20, method_list = c("mr_ivw"))
#' plot_mrcorge(res, scale='exp')
#' @export
plot_mrcorge <- function(mrcorge_result, scale = "linear") {
  plotdat <- mrcorge_result$result
  plotdat$category <- factor(plotdat$category, levels = c("Group-specific", "Cumulative"))
  plotdat$group <- factor(plotdat$group, levels = 1:max(plotdat$group))
  
  if (scale == "linear") {
    plotdat$b_plot <- plotdat$b
    plotdat$y_min <- plotdat$b - 1.96 * plotdat$se
    plotdat$y_max <- plotdat$b + 1.96 * plotdat$se
  } else if (scale == "exp") {
    plotdat$b_plot <- exp(plotdat$b)
    plotdat$y_min <- exp(plotdat$b - 1.96 * plotdat$se)
    plotdat$y_max <- exp(plotdat$b + 1.96 * plotdat$se)
  } else {
    stop("Please choose linear for continuous traits and exp for binary traits!")
  }
  
  base_plot <- ggplot(plotdat[plotdat$method != "MR Egger intercept", ], aes(x = group, y = b_plot, color = category)) +
    geom_hline(yintercept = ifelse(scale == "linear", 0, 1), lty = 2, col = "darkgrey") +
    facet_wrap(~method) +
    geom_point(position = position_dodge(width = 0.3)) +
    geom_errorbar(aes(ymin = y_min, 
                      ymax = y_max), width=0, position = position_dodge(width = 0.3)) +
    theme_minimal() +
    xlab("Instrument group (more likely core instruments => more likely peripheral instruments)") +
    ylab(ifelse(scale == "linear", expression(hat(beta)[MR]), expression(hat(OR)[MR]))) +
    labs(color = "") +
    scale_color_manual(values = c("orange", "darkblue")) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.position = "bottom")
  
  base_plot
}



