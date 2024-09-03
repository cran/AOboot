#' @title AOboot_one Class
#' @description A S3 class to represent one-way ANOVAs.
#' @param type.aov Character string giving the type of ANOVA computed.
#' @param factor Names of groups in the entered factor.
#' @param anova Results of the ANOVA.
#' @param em Bootstrapped estimated marginal means.
#' @param no.test Number of tests conducted that did not produce errors.
#' @param ph Bootstrapped post hoc tests.
#' @return An object of class \code{"AOboot_one"}.

AOboot_one <- function(type.aov,
                       factor,
                       anova,
                       em,
                       no.test,
                       ph){

  structure(list(type.aov,
                 factor,
                 anova,
                 em,
                 no.test,
                 ph),
            class = "AOboot_one")
}

#' @export print.AOboot_one
#' @export

print.AOboot_one <- function(x, ...){
  cat("\n", x$type.aov, "\n")
  cat("factor: ", toString(x$factor), "\n\n")
  print(x$anova)
  cat("\n\n",
          "***********************************************************",
          "\n\n")
  cat("Estimated marginal means \n")
  print(x$em)
  cat("\nPost-hoc tests\n")
  cat("number of tests: ", x$no.test, "\n")
  print(x$ph)
}
