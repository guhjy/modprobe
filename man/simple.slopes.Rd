\name{simple.slopes}
\alias{simple.slopes}
\alias{print.simple.slopes}
\alias{plot.simple.slopes}
\title{
Simple Slopes
}
\description{
\code{simple.slopes()} generates simple slopes and tests for the effect of a predictor on an outcome given a linear regression fit. \code{plot.simple.slopes()} generates a plot of the simple regression lines.
}
\usage{
simple.slopes(fit, predictor = NULL, moderator = NULL, at.mod.level = NULL, mod.level.names = NULL, alpha = .05, sig.region = NULL, ...)

\method{print}{simple.slopes}(s, ...)

\method{plot}{simple.slopes}(s, pred.range = NULL, colors = NULL, ...)
}
\arguments{
  \item{fit}{
the output of a call to \code{lm()} containing a linear interaction between a continuous predictor and a one or two moderators. If other variables are included in the model, all output will be taken with those variables set to 0.
}
  \item{predictor}{
optional; a string containing the name of the predictor. By default, it is assumed that the first variable to appear in fit is the predictor; if not, it must be specified here.
}
  \item{moderator}{
optional; a string containing the name of the moderator. By default, it is assumed that the second variable to appear in fit is the sole moderator; if not, it must be specified here. If more than one moderator is present, all of them must be named because the default is to only use the second variable in the fit.
}
  \item{at.mod.level}{
optional; a list of values of each moderator for which simple slopes are to be computed. If \code{NULL}, will compute simple slopes at the mean, mean - SD, and mean + SD of the moderator if the moderator is continuous or at the unique values of the moderator if it is categorical. Ideally, each item in the list will be named in accordance with each moderator, but otherwise each item will be used in the order it appears.
}
  \item{mod.level.names}{
optional; a vector of the same length as \code{at.mod.level} containing the names of the levels of the moderator for which simple slopes are to be calculated, in the same order as \code{at.mod.level}. If \code{NULL}, will use the values of the moderator as the names, unless \code{at.mod.level} is also \code{NULL} and the moderator is continuous, in which case "Mean", "Mean - SD", and "Mean + SD" will be used.
}
  \item{pred.range}{
optional; a vector of length 2 containing the range of the predictor for which the simple regression lines are to be plotted in a call to \code{plot.simple.slopes()}. If \code{NULL}, will use the range of the predictor given in \code{fit}.
}
  \item{alpha}{
the alpha value (significance criterion) of the t-tests for simple slopes.
}
  \item{sig.region}{
optional; a \code{sig.region} object. If specified, will compute and test the simple slopes at the significance bounds. This is generally not of use. If \code{at.mod.level} is specified, \code{sig.region} will be ignored.
}
  \item{s}{
a \code{simple.slopes} object; the output of a call to \code{simple.slopes()}.
}
  \item{colors}{
a vector of colors to be used to specify the colors of the lines in the simple slopes plot. See 'Color Specification' at \code{\link{par}}. 
}
  \item{...}{
other arguments passed to or from other functions.
}
}
\details{
\code{simple.slopes()} replicates the simple slopes tool at [quantpsy.org] by allowing users to enter a regression object and specify options for display and computation. Currently it can support up to two moderators. Categorical moderators are supported, but only if exactly one moderator is present. A categorical variable should be entered as a factor variable in the original fit.

It has its own print method, which produces the values and tests of the simple slopes and intercepts for the given values of the moderator(s). \code{plot()} generates a plot of the simple slopes. Note that this plot is a \code{ggplot} object, meaning users can edit it with \pkg{ggplot2} syntax.

}
\value{
\code{simple.slopes}: An object of class "simple.slopes" with the following entries:
\code{simple.lines} The slopes (w0) and intercepts (w1) and their t- and p-values and confidence intervals for each simple line.
\code{var.names} The names of the outcome, predictor, and moderator(s).
\code{pred.range} The range of the predictor in the original data set.
\code{cz} The values and names of the conditional values of the moderator(s) at which the simple slopes are evaluated.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Noah Greifer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#One Moderator
data("mpg")
fit <- lm(cty ~ displ * drv, data = mpg)
summary(fit)

s <- simple.slopes(fit)
s
plot(s)

#Two Moderators
fit2 <- lm()
}