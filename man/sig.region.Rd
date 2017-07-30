\name{sig.region}
\alias{sig.region}
\alias{print.sig.region}
\alias{plot.sig.region}
\title{
Johnson-Neyman Significance Regions
}
\description{
Estimates the Johnson-Neyman significance bounds for a linear fit. These are the bounds at which the effect of the predictor on the outcome changes from being significant to being non-signficant given varying values of the moderator(s). Also estmates the confidence bands for the slope of the predictor on the outcome for varying levels of the moderator(s).}
\usage{
sig.region(fit, predictor = NULL, moderator = NULL, mod.range = NULL, at.mod2.level = NULL, mod2.level.names = NULL, alpha = 0.05, simultaneous = FALSE)
}
\method{print}{sig.region}(sr, digits = 4, ...)

\method{plot}{sig.region}(sr, mod.range = NULL, colors = NULL, ggplot2 = TRUE, ...)
}\arguments{
  \item{fit}{
the output of a call to \code{lm()} containing a linear interaction between a continuous predictor and one or two moderators. If other variables are included in the model, all output will be taken with those variables set to 0.}
  \item{predictor}{
optional; a string containing the name of the predictor. By default, it is assumed that the first variable to appear in fit is the predictor; if not, it must be specified here.}
  \item{moderator}{
optional; a string containing the name of the moderator. By default, it is assumed that the second variable to appear in fit is the sole moderator; if not, it must be specified here. If more than one moderator is present, all of them must be named because the default is to only use the second variable in the fit.
}
  \item{mod.range}{
optional; a vector of length 2 containing the range of the first moderator for which the simple regression lines are to be plotted in a call to \code{plot.sig.region()}. The output of \code{sig.region()} will also note whether the computed significance bounds are within the range of the moderator. If \code{NULL}, will use the range of the moderator given in \code{fit}.
}
  \item{at.mod2.level}{
optional; a list of values of the second moderator (if present) for which the signfiicance bounds and regions are to be computed. If \code{NULL} and two moderators are specified, will compute the signfiicance bounds and regions at the mean, mean - SD, and mean + SD of the second moderator. Ideally, each item in the list will be named in accordance with each moderator, but otherwise each item will be used in the order it appears.
}
  \item{mod2.level.names}{
optional; a vector of the same length as \code{at.mod2.level} containing the names of the levels of the second moderator for which the signfiicance bounds and regions are to be calculated, in the same order as \code{at.mod2.level}. If \code{NULL}, will use the values of the second moderator as the names, unless \code{at.mod2.level} is also \code{NULL} and the moderator is continuous, in which case "Mean", "Mean - SD", and "Mean + SD" will be used.
}
  \item{alpha}{
the alpha value (significance criterion) corresponding to the confidence level for the signfiicance bounds and regions.
}
  \item{simultaneous}{
whether the bounds and bands should be calculated considering the tests as simultaneous or not. If so, a stricter criterion is used, making the bands wider.
}
  \item{sr}{
a \code{sig.region} object; the output of a call to \code{sig.region()}.
}
  \item{digits}{
the number of digits to be displayed for the significance bounds.
}
  \item{colors}{
a vector of colors to be used to specify the colors of the lines in the significance bands plot. See 'Color Specification' at \code{\link{par}}. 
}
  \item{ggplot2}{
logical; whether to print results using \pkg{ggplot2} or \R's \pkg{graphics}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}