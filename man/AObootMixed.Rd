\name{AObootMixed}
\alias{AObootMixed}
\title{Bootstrapped ANOVA for Mixed Designs}
\description{
In case of violations of the assumption of the normal distribution, researchers
usually employ bootstrapping. Based on the \code{R} packages
\href{https://cran.r-project.org/package=afex/}{afex} and
\href{https://cran.r-project.org/package=emmeans/}{emmeans}, this function
computes bootstrapped confidence intervals for the effect sizes, estimated
marginal means, and post hoc tests for Mixed ANOVAs. Furthermore, the p-values
of the F-statistic are adjusted to reflect the probability to obtain equal or
higher values than the raw, non-bootstrapped ANOVA (Stine, 1989 <doi:10.1177/0049124189018002003>; see
also \href{https://www.uvm.edu/~statdhtx/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html}{this tutorial by Nadine Spychala}.).}

\usage{
AObootMixed(var.within,
            var.between,
            var.id,
            levels.w1,
            levels.b1,
            eff.si = c("pes", "ges"),
            data,
            silence = FALSE,
            n.sim = 1000,
            alpha = .05,
            seed = 1234,
            n.round = 2)
            }

\arguments{
  \item{var.within}{Variable(s) reflecting the within-subject level.}
  \item{var.between}{Variable(s) reflecting the between-subject level.}
  \item{var.id}{Unique person specifier.}
  \item{levels.w1}{Levels of the within-subjects variable. Must be identical
  with the levels in the dataset.}
  \item{levels.b1}{Levels of the between-subjects variable. Must be identical
  with the levels in the dataset.}
  \item{eff.si}{Effect size for the F-tests. \code{pes} reflects partial
  eta-squared, \code{ges} reflects eta-squared.}
  \item{data}{Name of the dataframe. The dataset must be in a wide-format, with
  one row per participant.}
  \item{silence}{Logical. If FALSE, progress of the bootstrapping procedure will
  be displayed.}
  \item{n.sim}{Number of bootstrap samples to be drawn.}
  \item{alpha}{Type I error.}
  \item{seed}{To make the results reproducible, it is recommended to set a
  random seed parameter.}
  \item{n.round}{Number of digits in the output.}}

\details{The p-value of the F-test (\code{Pr(>F)}) in the output reflects the
probability to obtain an F-value as high as or higher than the F-value from the
raw, non-bootstrapped ANOVA. Thus, it should not be mistaken as a p-value in the
sense of a null hypothesis significance test. More information about this can be
found \href{https://www.uvm.edu/~statdhtx/StatPages/Resampling/Bootst1way/bootstrapping_oneway.html}{in this tutorial by Nadine Spychala}.}

\value{
\item{type.aov}{Type of ANOVA conducted.}
\item{factor1}{Name of the groups in the between factor.}
\item{factor2}{Name of the groups in the within factor.}
\item{anova}{Results of the conducted ANOVA (i.e., degrees of freedom, F-test,
p-value, effect size with bootstrap confidence interval, and numbers of tests
for which convergence was achieved.}
\item{em.1}{Estimated marginal means for between factor.}
\item{em.2}{Estimated marginal means for within factor.}
\item{em.3}{Estimated marginal means for between factor by within factor.}
\item{em.4}{Estimated marginal means for within factor by between factor.}
\item{no.test1}{Number of post hoc tests for the between factor for which
convergence was achieved.}
\item{no.test2}{Number of post hoc tests for the within factor for which
convergence was achieved.}
\item{no.test3}{Number of post hoc tests for the between factor by within factor
for which convergence was achieved.}
\item{no.test4}{Number of post hoc tests for within factor by between factor for
which convergence was achieved.}
\item{ph.1}{Post hoc tests for between factor.}
\item{ph.2}{Post hoc tests for within factor.}
\item{ph.3}{Post hoc tests for between factor by within factor.}
\item{ph.4}{Post hoc tests for within factor by between factor.}


  output <- list(type.aov = "Two-way mixed ANOVA",
                 factor1 = levels.b1,
                 factor2 = levels.w1,
                 anova = round(orig.aov$anova_table, n.round),
                 em.1 = dat.em1,
                 no.test1 = no.test1,
                 ph.1 = dat.ph1,
                 em.2 = dat.em2,
                 no.test2 = no.test2,
                 ph.2 = dat.ph2,
                 em.3 = dat.em3,
                 no.test3 = no.test3,
                 ph.3 = dat.ph3,
                 em.4 = dat.em4,
                 no.test4 = no.test4,
                 ph.4 = dat.ph4)

}

\references{
Stine, R. (1989). An introduction to bootstrap methods: Examples and ideas.
Sociological Methods & Research, 18(2-3), 243--291. <https://doi.org/10.1177/0049124189018002003>}

\examples{
library(carData)

# The OBrienKaiser dataset from the carData package

ao <- OBrienKaiser

# Add a unique person identifier to the dataset

ao$pers <- 1:nrow(OBrienKaiser)

# Mixed ANOVA
\donttest{
AObootMixed(
    var.within = c("pre.1", "post.1", "fup.1"),
    var.between = "treatment",
    var.id = "pers",
    levels.w1 = c("pre", "post", "fup"),
    levels.b1 = c("control", "A", "B"),
    eff.si = "pes",
    data = ao,
    n.sim = 1000,
    alpha = .05,
    seed = 1234,
    n.round = 2)
}
}

\author{
Lisa-Marie Segbert, Christian Blötner
\email{c.bloetner@gmail.com}}

\keyword{methods}
\keyword{htest}
\keyword{nonparametric}
