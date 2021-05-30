
anova_summary <- function(object, effect.size = "ges",  observed = NULL){
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
  .args <- attr(object, "args") # exist only in anova_test()
  results <- add_anova_effect_size(results)

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
  aov.table <- .summary$univariate.tests %>%
    convert_anova_object_as_data_frame() %>%
    set_colnames(c("Effect", "SSn", "DFn", "SSd", "DFd", "F", "p")) %>%
    select(
      .data$Effect, .data$DFn, .data$DFd,
      .data$SSn, .data$SSd, .data$F, .data$p
    ) %>%
    mutate(`p<.05` = ifelse(.data$p < 0.05, "*",''))
  sphericity.test <- corrections <- NULL
  # Mauchly's Test for Sphericity
  if(nrow(.summary$sphericity.tests) > 0){
    sphericity.test <- .summary$sphericity.tests %>%
      convert_anova_object_as_data_frame() %>%
      set_colnames(c("Effect", "W", "p")) %>%
      mutate(`p<.05` = ifelse(.data$p < 0.05, "*",''))
  }
  # Sphericity corrections
  if(nrow(.summary$sphericity.tests) > 0){
    corrections <- .summary$pval.adjustments %>%
      as.data.frame() %>%
      set_colnames(c("GGe", "p[GG]", "HFe", "p[HF]")) %>%
      tibble::rownames_to_column("Effect")
    p.gg.signif <- ifelse(corrections[["p[GG]"]] < 0.05, "*",'')
    p.hf.signif <- ifelse(corrections[["p[HF]"]] < 0.05, "*",'')
    corrections <- corrections %>%
      add_column(`p[GG]<.05` = p.gg.signif, .after = "p[GG]") %>%
      add_column(`p[HF]<.05` = p.hf.signif, .after = "p[HF]")
  }
  # Results
  results <- list(ANOVA = aov.table)
  if(!is.null(sphericity.test)){
    results $`Mauchly's Test for Sphericity` <- sphericity.test
    results$`Sphericity Corrections` <-  corrections
    results <- results %>% add_corrected_df()
  }
  results
}

convert_anova_object_as_data_frame <- function(aov.table){
  aov.table.list <- list(Effect = rownames(aov.table))
  for(col in colnames(aov.table)){
    aov.table.list[[col]] <- aov.table[, col]
  }
  aov.table <- as.data.frame(aov.table.list, stringsAsFactors = FALSE)
  rownames(aov.table) <- 1:nrow(aov.table)
  aov.table
}

add_corrected_df <- function(.summary){
  aov.table <- .summary$ANOVA %>%
    select(.data$Effect, .data$DFn, .data$DFd)
  corrections <- .summary$`Sphericity Corrections` %>%
    dplyr::left_join(aov.table, by = "Effect") %>%
    mutate(
      df.gg = paste(round_value(.data$GGe*.data$DFn, 2), round_value(.data$GGe*.data$DFd, 2), sep = ", "),
      df.hf = paste(round_value(.data$HFe*.data$DFn, 2), round_value(.data$HFe*.data$DFd, 2), sep = ", ")
    ) %>%
    select(-.data$DFd, -.data$DFn)
  df.gg <- corrections$df.gg
  df.hf <- corrections$df.hf
  .summary$`Sphericity Corrections` <- corrections %>%
    select(-.data$df.gg, -.data$df.hf) %>%
    add_column(`DF[GG]` = df.gg, .after = "GGe") %>%
    add_column(`DF[HF]` = df.hf, .after = "HFe")
  .summary
}

# Summary of independent anova
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_independent_anova <- function(res.anova){
  res.anova <- as.data.frame(res.anova)
  .residuals <- res.anova["Residuals", 1:2]
  if('Mean Sq' %in% colnames(res.anova)){
    # exists when res.anova is from stats::anova
    res.anova <- select(res.anova, -.data$`Mean Sq`)
  }
  if('Sum Sq' %in% colnames(res.anova)){
    # in stats::anova, Sum Sq is not the first column, so do select
    res.anova <- res.anova %>% select(.data$`Sum Sq`, dplyr::everything())
    colnames(res.anova) <- c('SSn','DFn','F','p')
    ss.exists <- TRUE
  }
  else{
    # case of white.adjust = TRUE. SS doesnt exist in the results
    colnames(res.anova) <- c('DFn','F','p')
    ss.exists <- FALSE
  }
  res.anova <- res.anova %>%
    tibble::rownames_to_column("Effect")  %>%
    add_column(DFd = .residuals$Df, .after = "DFn") %>%
    mutate(`p<.05` = ifelse(.data$p < 0.05, "*",'')) %>%
    filter(.data$Effect != "Residuals")
  if(ss.exists){
    res.anova <- res.anova %>%
      add_column(SSd = .residuals$`Sum Sq`, .after = "SSn")
  }
  results <- list(ANOVA = res.anova)
  results
}

# Summary of anova from stats::aov
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
summary_aov <- function(res.anova){
  remove_empty_space <- function(x){
    sapply(x, function(x){strsplit(x, " ")[[1]][1]})
  }
  reformat_aov_summary <- function(aov.summary){
    if(inherits(aov.summary, "listof")){
      aov.summary <- as.data.frame(aov.summary[[1]])
    } else {as.data.frame(aov.summary)}
    .residuals <- aov.summary["Residuals", 1:2]
    aov.summary <- aov.summary %>%
      set_colnames(c("DFn", "SSn", "MS", "F", "p")) %>%
      tibble::rownames_to_column("Effect")  %>%
      add_column(DFd = .residuals$Df, .after = "DFn") %>%
      add_column(SSd = .residuals$`Sum Sq`, .after = "SSn") %>%
      mutate(`p<.05` = as.character(ifelse(.data$p < 0.05, "*",''))) %>%
      mutate(Effect = remove_empty_space(.data$Effect)) %>%
      filter(!is.na(.data$p)) %>%
      select(-.data$MS)
    aov.summary
  }
  res.anova <- summary(res.anova) %>%
    map(reformat_aov_summary) %>%
    dplyr::bind_rows() %>%
    order_by_interaction_levels()
  results <- list(ANOVA = res.anova)
  results
}


# Reorder ANOVA table by interaction levels in the term
#order_by_interaction_levels <- function(aov.table){
#  .terms <- aov.table$Effect
#  nb.interaction <- str_count(.terms, ":")
#  aov.table %>% dplyr::arrange(nb.interaction)
#}


# Add effect size
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
add_anova_effect_size <- function(res.anova.summary){
  ss.exists <- "SSn" %in% colnames(res.anova.summary$ANOVA)
  if (!ss.exists) {
    return(res.anova.summary)
  }
  res.anova.summary <- add_partial_eta_squared(res.anova.summary)
  return(res.anova.summary)
}


# Partial eta squared
add_partial_eta_squared <- function(res.anova.summary){
  res.anova.summary$ANOVA$`.data`$pes <- res.anova.summary$ANOVA$`.data`$SSn/(res.anova.summary$ANOVA$`.data`$SSn + res.anova.summary$ANOVA$`.data`$SSd) #%>%
    #mutate(pes = .data$SSn/(.data$SSn + .data$SSd))
  return(res.anova.summary)
}



