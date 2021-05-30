
anova_summary <- function(object){
  if(inherits(object, "Anova.mlm")){
    results <- repeated_anova_summary(object)
  }
  else if(inherits(object, "anova")){
    results <- summary_independent_anova(object)
  }
  else if(inherits(object, c("aov", "aovlist"))){
    results <- summary_aov(object)
  }
  else{
    stop("Non-supported object passed: ",
         paste(class(object), collapse = ", "), ". ",
         "Object needs to be of class 'Anova.mlm' or 'anova'.")
  }
  #.args <- attr(object, "args") # exist only in anova_test()
  #results <- add_anova_effect_size(results)

  results$pes = results$F.value * results$df1 / (results$F.value * results$df1 + results$df2)

  #if(!detailed){
  #  results <- remove_details(results, method = "anova")
  #}
  #results$ANOVA <- order_by_interaction_levels(results$ANOVA)
  #results <- results %>% map(~dplyr::mutate_if(., is.numeric, round_value, 3))
  #if(length(results) == 1) results <- results[[1]]
  #results %>% set_attrs(args = .args)
  return(results)
}




# Summary of Anova.mlm object: summary_anova_mlm
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# this function is used for repeated and mixed anova
repeated_anova_summary <- function(res.anova){
  .summary <- suppressWarnings(summary(res.anova))
  # Anova table converted into data frame
  #aov.table <- .summary$univariate.tests %>%
  #  convert_anova_object_as_data_frame() %>%
  #  set_colnames(c("Effect", "SSn", "DFn", "SSd", "DFd", "F", "p")) %>%
  #  select(
  #    .data$Effect, .data$DFn, .data$DFd,
  #    .data$SSn, .data$SSd, .data$F, .data$p
  #  ) %>%
  #  mutate(`p<.05` = ifelse(.data$p < 0.05, "*",''))

  aov.table = convert_anova_object_as_data_frame(.summary$univariate.tests)

  colnames(aov.table) = c("Effect", "SSn", "df1", "SSd", "df2", "F.value", "p.value")

  aov.table = aov.table[c("Effect","df1", "df2","SSn","SSd","F.value","p.value")]

  #aov.table$`p<.05` = ifelse(aov.table$p < 0.05, "*",'')

  #sphericity.test <- corrections <- NULL
  # Mauchly's Test for Sphericity
  #if(nrow(.summary$sphericity.tests) > 0){
  #  #sphericity.test <- .summary$sphericity.tests %>%
  #  #  convert_anova_object_as_data_frame() %>%
  #  #  set_colnames(c("Effect", "W", "p")) %>%
  #  #  mutate(`p<.05` = ifelse(.data$p < 0.05, "*",''))
  #  sphericity.test <- convert_anova_object_as_data_frame(.summary$sphericity.tests)
  #  colnames(sphericity.test) = c("Effect", "W", "p")
  #  sphericity.test$`p<.05` = ifelse(sphericity.test$p < 0.05, "*",'')
  #}
  # Sphericity corrections
  #if(nrow(.summary$sphericity.tests) > 0){
  #  #corrections <- .summary$pval.adjustments %>%
  #  #  as.data.frame() %>%
  #  #  set_colnames(c("GGe", "p[GG]", "HFe", "p[HF]")) %>%
  #  #  tibble::rownames_to_column("Effect")
  #  corrections <- as.data.frame(.summary$pval.adjustments)
  #  colnames(corrections) = c("GGe", "p[GG]", "HFe", "p[HF]")
  #  corrections$Effect = row.names(corrections)
  #  row.names(corrections) = 1:nrow(corrections)
  #  p.gg.signif <- ifelse(corrections[["p[GG]"]] < 0.05, "*",'')
  #  p.hf.signif <- ifelse(corrections[["p[HF]"]] < 0.05, "*",'')
  #
  # #corrections <- corrections %>%
  #  #  add_column(`p[GG]<.05` = p.gg.signif, .after = "p[GG]") %>%
  #  #  add_column(`p[HF]<.05` = p.hf.signif, .after = "p[HF]")
  #
  #  corrections$`p[GG]<.05` = p.gg.signif
  #  corrections$`p[HF]<.05` = p.hf.signif
  #}
  # Results
  #results <- list(ANOVA = aov.table)
  #if(!is.null(sphericity.test)){
  #  results $`Mauchly's Test for Sphericity` <- sphericity.test
  #  results$`Sphericity Corrections` <-  corrections
  #  results <- add_corrected_df(results)
  #}
  aov.table = aov.table[c("Effect", "df1", "df2", "F.value", "p.value")]
  return(aov.table)
}

convert_anova_object_as_data_frame <- function(aov.table){
  aov.table.list <- list(Effect = rownames(aov.table))
  for(col in colnames(aov.table)){
    aov.table.list[[col]] <- aov.table[, col]
  }
  aov.table <- as.data.frame(aov.table.list, stringsAsFactors = FALSE)
  rownames(aov.table) <- 1:nrow(aov.table)
  return(aov.table)
}

# Summary of independent anova
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_independent_anova <- function(res.anova){
  res.anova <- as.data.frame(res.anova)
  .residuals <- res.anova["Residuals", 1:2]
  if('Mean Sq' %in% colnames(res.anova)){
    # exists when res.anova is from stats::anova
    # res.anova <- select(res.anova, -.data$`Mean Sq`)
    dropvar1 = names(res.anova) %in% c("Mean Sq")
    res.anova = res.anova[!dropvar1]
  }
  if('Sum Sq' %in% colnames(res.anova)){
    # in stats::anova, Sum Sq is not the first column, so do select
    #res.anova <- res.anova %>% select(.data$`Sum Sq`, dplyr::everything())
    res.anova = res.anova[c("Sum Sq","Df","F value","Pr(>F)")]
    colnames(res.anova) <- c('SSn','df1','F','p')
    ss.exists <- TRUE
  }else{
    # case of white.adjust = TRUE. SS doesnt exist in the results
    colnames(res.anova) <- c('df1','F','p')
    ss.exists <- FALSE
  }
  res.anova$Effect = row.names(res.anova)
  row.names(res.anova) = 1:nrow(res.anova)
  res.anova$df2 = .residuals$Df
  #res.anova$`p<.05` = ifelse(res.anova$p < 0.05,
  #                           "*",'')
  res.anova = subset(res.anova,
                     Effect != "Residuals")
  #res.anova <- res.anova %>%
  #  tibble::rownames_to_column("Effect")  %>%
  #  add_column(DFd = .residuals$Df, .after = "DFn") %>%
  #  mutate(`p<.05` = ifelse(.data$p < 0.05, "*",'')) %>%
  #  filter(.data$Effect != "Residuals")
  #if(ss.exists){
  #  res.anova$SSd <- .residuals$`Sum Sq`
  #    #add_column(SSd = .residuals$`Sum Sq`, .after = "SSn")
  #}
  res.anova = res.anova[c("Effect", "df1", "df2", "F", "p")]
  colnames(res.anova) = c("Effect", "df1", "df2", "F.value", "p.value")
  return(res.anova)
  #results <- list(ANOVA = res.anova)
  #results
}

# Summary of anova from stats::aov
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_aov <- function(res.anova){
  remove_empty_space <- function(x){
    sapply(x, function(x){strsplit(x, " ")[[1]][1]})
  }

  #res.anova <- summary(res.anova) %>%
  #  map(reformat_aov_summary) %>%
  #  dplyr::bind_rows() %>%
  #  order_by_interaction_levels()
  #results <- list(ANOVA = res.anova)
  #results

  res.anova = summary(res.anova)
  if(length(res.anova) != 1){
    res.anova = lapply(res.anova, reformat_aov_summary)
    res.anova = do.call(rbind, res.anova)
  } else{
    res.anova = reformat_aov_summary(res.anova)
  }
  res.anova = subset(res.anova,
                     !is.na(p.value))
  row.names(res.anova) = 1:nrow(res.anova)
  res.anova = res.anova[c("Effect", "df1", "df2", "F.value", "p.value")]
  return(res.anova)
}

reformat_aov_summary <- function(aov.summary){
  if(inherits(aov.summary, "listof")){
    aov.summary <- as.data.frame(aov.summary[[1]])
  } else {as.data.frame(aov.summary)}
  .residuals <- aov.summary["Residuals", 1:2]
  #aov.summary <- aov.summary %>%
  #  set_colnames(c("df1", "SSn", "MS", "F.value", "p.value")) %>%
  #  tibble::rownames_to_column("Effect")  %>%
  #  add_column(df2 = .residuals$Df, .after = "df1") %>%
  #  add_column(SSd = .residuals$`Sum Sq`, .after = "SSn") %>%
  #  mutate(`p<.05` = as.character(ifelse(.data$p < 0.05, "*",''))) %>%
  #  mutate(Effect = remove_empty_space(.data$Effect)) %>%
  #  filter(!is.na(.data$p)) %>%
  #  select(-.data$MS)

  colnames(aov.summary) = c("df1", "SSn", "MS", "F.value", "p.value")
  aov.summary$Effect = row.names(aov.summary)
  row.names(aov.summary) = 1:nrow(aov.summary)
  aov.summary$df2 = .residuals$Df
  aov.summary = aov.summary[c("Effect", "df1", "df2", "F.value", "p.value")]
  return(aov.summary)
}


# Reorder ANOVA table by interaction levels in the term
#order_by_interaction_levels <- function(aov.table){
#  .terms <- aov.table$Effect
#  nb.interaction <- str_count(.terms, ":")
#  aov.table %>% dplyr::arrange(nb.interaction)
#}


# Add effect size
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#add_anova_effect_size <- function(res.anova.summary){
#  ss.exists <- "SSn" %in% colnames(res.anova.summary$ANOVA)
#  if (!ss.exists) {
#    return(res.anova.summary)
#  }
#  res.anova.summary <- add_partial_eta_squared(res.anova.summary)
#  return(res.anova.summary)
#}


# Partial eta squared
#add_partial_eta_squared <- function(res.anova.summary){
#  res.anova.summary$ANOVA$`.data`$pes <- res.anova.summary$ANOVA$`.data`$SSn/(res.anova.summary$ANOVA$`.data`$SSn + res.anova.summary$ANOVA$`.data`$SSd) #%>%
#    #mutate(pes = .data$SSn/(.data$SSn + .data$SSd))
#  return(res.anova.summary)
#}



