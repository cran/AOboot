#' @title AOboot_two Class
#' @description A S3 class to represent two-way ANOVAs.
#' @param type.aov Character string giving the type of ANOVA computed.
#' @param factor1 Names of groups in the first factor.
#' @param factor2 Names of groups in the second factor.
#' @param anova Results of the ANOVA.
#' @param em.1 Bootstrapped estimated marginal means for factor 1.
#' @param no.test1 Number of bootstrapped tests conducted for factor 1 that did
#' not produce errors.
#' @param ph.1 Bootstrapped post hoc tests for factor 1.
#' @param em.2 Bootstrapped estimated marginal means for factor 2.
#' @param no.test2 Number of bootstrapped tests conducted for factor 2 that did
#' not produce errors.
#' @param ph.2 Bootstrapped post hoc tests for factor 2.
#' @param em.3 Bootstrapped estimated marginal means for factor 1 by factor 2.
#' @param no.test3 Number of bootstrapped tests conducted for factor 1 by factor
#' 2 that did not produce errors.
#' @param ph.3 Bootstrapped post hoc tests for factor 1 by factor 2.
#' @param em.4 Bootstrapped estimated marginal means for factor 2 by factor 1.
#' @param no.test4 Number of bootstrapped tests conducted for factor 2 by factor
#' 1 that did not produce errors.
#' @param ph.4 Bootstrapped post hoc tests for factor 2 by factor 1.
#' @return An object of class \code{"AOboot_two"}.

AOboot_two <- function(type.aov,
                       factor1,
                       factor2,
                       anova,
                       em.1,
                       no.test1,
                       ph.1,
                       em.2,
                       no.test2,
                       ph.2,
                       em.3,
                       no.test3,
                       ph.3,
                       em.4,
                       no.test4,
                       ph.4){

  structure(list(type.aov,
                 factor1,
                 factor2,
                 anova,
                 em.1,
                 no.test1,
                 ph.1,
                 em.2,
                 no.test2,
                 ph.2,
                 em.3,
                 no.test3,
                 ph.3,
                 em.4,
                 no.test4,
                 ph.4),
            class = "AOboot_two")
}

#' @export print.AOboot_two
#' @export

print.AOboot_two <- function(x, ...){
  cat("\n", x$type.aov, "\n")
  cat("factor 1: ", toString(x$factor1), "\n")
  cat("factor 2: ", toString(x$factor2), "\n\n")
  print(x$anova)
  cat("\n\n",
          "*************************************************************",
          "\n\n")
  cat("Estimated marginal means factor 1\n")
  print(x$em.1)
  cat("\nPost-hoc tests factor 1\n")
  cat("number of tests: ", x$no.test1, "\n")
  print(x$ph.1)
  cat("\n\n",
          "*************************************************************",
          "\n\n")
  cat("Estimated marginal means factor 2\n")
  print(x$em.2)
  cat("\nPost-hoc tests factor 2\n")
  cat("number of tests: ", x$no.test2, "\n")
  print(x$ph.2)
  cat("\n\n",
          "*************************************************************",
          "\n\n")
  cat("Estimated marginal means factor 1 by factor 2\n")
  print(x$em.3)
  cat("\nPost-hoc tests factor 1 by factor 2\n")
  cat("number of tests: ", x$no.test3, "\n")
  print(x$ph.3)
  cat("\n\n",
          "*************************************************************",
          "\n\n")
  cat("Estimated marginal means factor 2 by factor 1\n")
  print(x$em.4)
  cat("\nPost-hoc tests factor 2 by factor 1\n")
  cat("number of tests: ", x$no.test4, "\n")
  print(x$ph.4)
}
