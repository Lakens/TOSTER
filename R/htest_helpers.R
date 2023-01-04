#' Helpers for \code{htest} objects
#'
#' Functions to help interpret or display objects of the class \code{htest}
#'
#' @param htest A S3 object of the class \code{htest}
#' @param test_statistics A logical variable to display the test statistics.
#' @param show_ci A logical variable to display the confidence interval.
#' @param extract_names A logical variable to take the names from the S3 object (i.e., statistic for \code{t.test} would be "t")
#'
#' @examples
#' # To be added
#'
#' @name htest-helpers

# Possible functions
# t.test wilcox.test oneway.test kruskal.test friedman.test
# cor.test

#' @rdname htest-helpers
#' @export
df_htest = function(htest,
                    test_statistics = TRUE,
                    show_ci = TRUE,
                    extract_names = TRUE){

  if(!("htest" %in% class(htest))){
    stop("htest must be of the class htest")
  }

  df1 = data.frame(method = htest$method)

  if(test_statistics) {


  # Get columns
  if(!is.null(htest$statistic)){
    if(length(htest$statistic) == 1){
      df1[ifelse(extract_names, names(htest$statistic),"statistic")] <- unname(htest$statistic)
    } else{
      for(i in 1:length(htest$statistic)){
        parm_name = ifelse(extract_names, names(htest$statistic)[i],paste0("statistic",i))
        df1[parm_name] <- unname(htest$statistic)[i]
      }
    }
  }

  if(!is.null(htest$parameter)){
    if(length(htest$parameter) == 1){
      df1[ifelse(extract_names, names(htest$parameter),"parameter")] <- unname(htest$parameter)
    } else{
      for(i in 1:length(htest$parameter)){
        parm_name = ifelse(extract_names, names(htest$parameter)[i],paste0("parameter",i))
        df1[parm_name] <- unname(htest$parameter)[i]
      }
    }
  }

  if(!is.null(htest$p.value)){
    df1$p.value <- unname(htest$p.value)
  }

  }

  if(!is.null(htest$estimate)){
    if(grepl("two sample t-test",htest$method,ignore.case = TRUE) && length(htest$estimate) > 1){
      htest$estimate = htest$estimate[1] - htest$estimate[2]
      names(htest$estimate) = c("mean difference")
    }
    if(length(htest$estimate) == 1){
      df1[ifelse(extract_names, names(htest$estimate),"estimate")] <- unname(htest$estimate)
    } else{
      for(i in 1:length(htest$estimate)){
        parm_name = ifelse(extract_names, names(htest$estimate)[i],paste0("estimate",i))
        df1[parm_name] <- unname(htest$estimate)[i]
      }
    }
  }

  if(!is.null(htest$stderr)){
    if(length(htest$stderr) == 1){
      df1["SE"] <- unname(htest$stderr)
    } else{
      for(i in 1:length(htest$stderr)){
        parm_name = paste0("SE",i)
        df1[parm_name] <- unname(htest$stderr)[i]
      }
    }
  }
  if(show_ci == TRUE){

  if(!is.null(htest$conf.int)){

    df1$lower.ci = min(htest$conf.int)
    df1$upper.ci = max(htest$conf.int)
    df1$conf.level = attr(htest$conf.int, "conf.level")

  }

}
  if(!is.null(htest$alternative)){
    df1$alternative = htest$alternative
  }

  if(!is.null(htest$null.value)){

    if(length(htest$null.value) == 1){
      df1["null"] <- unname(htest$null.value)
    } else{
      for(i in 1:length(htest$null.value)){
        parm_name = paste0("null",i)
        df1[parm_name] <- unname(htest$null.value)[i]
      }
    }
  }
  # Rename
  # names(df1)[names(df1) == 'old.var.name'] <- 'new.var.name'

  df_final = df1
  return(df_final)

}

