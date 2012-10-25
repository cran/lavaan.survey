\name{subset.survey.design}
\alias{subset.survey.design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a svydesign object consisting of a particular subset.
}
\description{
	For internal use only.
}
\usage{
subset.survey.design(survey.design, subset.call)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{survey.design}{
%%     ~~Describe \code{survey.design} here~~
}
  \item{subset.call}{
%%     ~~Describe \code{subset.call} here~~
}
}
\author{
	Daniel Oberski <daniel.oberski@gmail.com>
}

\seealso{
	\code{\link{lavaan.survey}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (survey.design, subset.call) 
{
    log <- eval(substitute(subset.call), envir = survey.design$variables, 
        enclos = parent.frame())
    survey.design$cluster <- subset(survey.design$cluster, log)
    survey.design$strata <- subset(survey.design$strata, log)
    survey.design$prob <- survey.design$prob[log]
    survey.design$allprob <- subset(survey.design$allprob, log)
    survey.design$variables <- subset(survey.design$variables, 
        log)
    survey.design$fpc$sampsize <- subset(survey.design$fpc$sampsize, 
        log)
    survey.design
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
