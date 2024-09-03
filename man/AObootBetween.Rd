\name{AObootBetween}
\alias{AObootBetween}
\title{Bootstrapped ANOVA for Between-Subject Designs}
\description{
In case of violations of the assumption of the normal distribution, researchers
usually employ bootstrapping. Based on the \code{R} packages
\href{https://cran.r-project.org/package=afex/}{afex} and
\href{https://cran.r-project.org/package=emmeans/}{emmeans}, this
function computes bootstrapped confidence intervals for the effect sizes,
estimated marginal means, and post hoc tests for one-way and two-way ANOVAs
following a between-subject design. Furthermore, the p-values of the F-statistic
are adjusted to reflect the probability to obtain equal or higher values than
the raw, non-bootstrapped ANOVA (Stine, 1989 <doi:10.1177/0049124189018002003>;
see also \href{https://www.uvm.edu/~statdhtx/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html}{this tutorial by Nadine Spychala}.).}

\usage{
AObootBetween(var.between,
              var.dv,
              var.id,
              levels.b1,
              levels.b2 = NULL,
              eff.si = c("pes", "ges"),
              data,
              silence = FALSE,
              n.sim = 1000,
              alpha = .05,
              seed = 1234,
              n.round = 2)
              }
\arguments{
  \item{var.between}{Variable(s) reflecting the between-subject level.}
  \item{var.dv}{Dependent variable.}
  \item{var.id}{Unique person specifier.}
  \item{levels.b1}{Levels of the first-named independent variable. Must be
  identical with the levels in the dataset.}
  \item{levels.b2}{For two-way ANOVAs. Levels of the second-named independent
  variable. Must be identical with the levels in the dataset.}
  \item{eff.si}{Effect size for the F-tests. "pes" reflects partial eta-squared,
  "ges" reflects eta-squared.}
  \item{data}{Name of the dataframe.}
  \item{silence}{Logical. If FALSE, progress of the bootstrapping procedure will
  be displayed.}
  \item{n.sim}{Number of bootstrap samples to be drawn.}
  \item{alpha}{Type I error.}
  \item{seed}{To make the results reproducible, it is recommended to set a
  random seed parameter.}
  \item{n.round}{Number of digits in the output.}}

\details{The p-value of the F-test (`Pr(>F)`) in the output reflects the
probability to obtain an F-value as high as or higher than the F-value from the
raw, non-bootstrapped ANOVA. Thus, it should not be mistaken as a p-value in the
sense of a null hypothesis significance test. More information about this can be
found \href{https://www.uvm.edu/~statdhtx/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html}{in this tutorial by Nadine Spychala}.}

\value{
\item{type.aov}{Type of ANOVA conducted.}
\item{factor}{Name of the groups in the factor (in one-way ANOVA).}
\item{factor1}{Name of the groups in the first factor (in two-way ANOVA).}
\item{factor2}{Name of the groups in the second factor (in two-way ANOVA).}
\item{anova}{Results of the conducted ANOVA (i.e., degrees of freedom, F-test,
p-value, effect size with bootstrap confidence interval, and numbers of tests
for which convergence was achieved.}
\item{em}{Estimated marginal means in one-way ANOVA.}
\item{em.1}{Estimated marginal means for factor 1 in two-way ANOVA.}
\item{em.2}{Estimated marginal means for factor 2 in two-way ANOVA.}
\item{em.3}{Estimated marginal means for factor 1 by factor 2 in two-way ANOVA.}
\item{em.4}{Estimated marginal means for factor 2 by factor 1 in two-way ANOVA.}
\item{no.test}{Number of post hoc tests in one-way ANOVAs for which convergence
was achieved.}
\item{no.test1}{Number of post hoc tests for factor 1 in two-way ANOVAs for
which convergence was achieved.}
\item{no.test2}{Number of post hoc tests for factor 2 in two-way ANOVAs for
which convergence was achieved.}
\item{no.test3}{Number of post hoc tests for factor 1 by factor 2 in two-way
ANOVAs for which convergence was achieved.}
\item{no.test4}{Number of post hoc tests for factor 2 by factor 1 in two-way
ANOVAs for which convergence was achieved.}
\item{ph}{Post hoc tests in one-way ANOVAs.}
\item{ph.1}{Post hoc tests for factor 1 in two-way ANOVAs.}
\item{ph.2}{Post hoc tests for factor 2 in two-way ANOVAs.}
\item{ph.3}{Post hoc tests for factor 1 by factor 2 in two-way ANOVAs.}
\item{ph.4}{Post hoc tests for factor 2 by factor 1 in two-way ANOVAs.}
}

\references{
Stine, R. (1989). An introduction to bootstrap methods: Examples and ideas.
Sociological Methods & Research, 18(2-3), 243--291. <https://doi.org/10.1177/0049124189018002003>}

\author{
Lisa-Marie Segbert, Christian Blötner
\email{c.bloetner@gmail.com}}

\examples{
library(carData)

# The OBrienKaiser dataset from the carData package

ao <- OBrienKaiser

# Add a unique person identifier to the dataset

ao$pers <- 1:nrow(OBrienKaiser)

# One-way between-subjects ANOVA
\donttest{
AObootBetween(
  var.between = "treatment",
  var.dv = "pre.1",
  var.id = "pers",
  levels.b1 = c("control", "A", "B"),
  eff.si = "ges",
  data = ao,
  n.sim = 1000,
  alpha = .05,
  seed = 1234,
  n.round = 2)
  }

# Two-way between-subjects ANOVA
\donttest{
AObootBetween(
  var.between = c("treatment", "gender"),
  var.dv = "pre.1",
  var.id = "pers",
  levels.b1 = c("control", "A", "B"),
  levels.b2 = c("M", "F"),
  eff.si = "pes",
  data = ao,
  n.sim = 1000,
  alpha = .05,
  seed = 1236,
  n.round = 2)
}
}

\keyword{methods}
\keyword{htest}
\keyword{nonparametric}
