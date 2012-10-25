\name{lavaan.survey}
\alias{lavaan.survey}
\title{
	Complex survey analysis of structural equation models (SEM)
}
\description{
	Takes a lavaan fit object and a complex survey design object as input
	and returns a structural equation modeling analysis based on the fit 
	object, where the complex sampling design is taken into account. 
	
	The structural equation model parameter estimates and standard errors
	are design-based. See Satorra and Muthen (1995) for details on the 
	procedure.
}
\usage{
lavaan.survey(lavaan.fit, survey.design)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lavaan.fit}{
	A \code{\linkS4class{lavaan}} object resulting from a lavaan call. 
	
	It probably makes most sense to use estimator="MLM" in the call
	since this is the estimator that will be used in the complex sample
	estimates, but it is not strictly necessary.
}
  \item{survey.design}{
	An  \code{\link{svydesign}} object resulting from a call to 
	svydesign in the survey package. This allows for incorporation of
	clustering, stratification, unequal probability weights, and finite
	population correction. See the survey documentation for more 
	information.
}
}
\details{
	The user specifies a complex sampling design with the survey package's
	\code{\link{svydesign}} function, and a structural equation model with
	\code{\link{lavaan}}.
	
	When calling lavaan.survey, the following steps are then completed:
	\enumerate{
		\item The covariance matrix of the observed variables
		(or matrices in the case of multiple
		group analysis) is estimated using the svyvar command from the
		survey package. 
	       	\item The asymptotic covariance matrix of the variances and
		covariances is obtained from the svyvar output  (the "Gamma"
		matrix)
		\item The lavaan model is re-fit using Maximum Likelihood
		with the covariance matrix as data. After normal-theory ML
		estimation, the standard errors (vcov matrix), likelihood ratio
		("chi-square") statistic, and all derived fit indices and
		statistics are adjusted for the complex sampling design using
		the Gamma matrix. I.e. the Satorra-Bentler (SB) corrections are
		obtained ("MLM" estimation in lavaan terminology).
	}

	The Satorra-Bentler ("aggregrated modeling") approach to complex 
	survey analysis of SEM was discussed by Satorra and Muthen (1995). 
	An alternative method to take clustering into account is multilevel 
	SEM ("disaggregated modeling").
}
\value{
	An object of class \code{\linkS4class{lavaan}}, where the estimates, 
	standard errors, vcov matrix, chi-square statistic, and fit measures 
	based on the chi-square take into account the complex survey 
	design. Several methods are available for \code{\linkS4class{lavaan}} 
	objects, including a \code{summary} method.}
\references{
	Oberski, D. and Saris, W. (2012). A model-based procedure to evaluate
	    the relative effects of different TSE components on structural equation
	    model parameter estimates. Presentation given at the International
	    Total Survey Error Workshop in Santpoort, the Netherlands. 
	    \url{http://daob.org/}

	Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics
		and standard errors in covariance structure analysis. 

	Satorra, A., and Muthen, B. (1995). Complex sample data in structural
	   equation modeling. Sociological methodology, 25, 267-316.
}
\author{
	Daniel Oberski - \url{http://daob.org} - \email{daniel.oberski@gmail.com}
}
\note{
    The function has been testing using simulation.
    
    Currently only continuous observed variables are implemented.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{svydesign}}
	\code{\link{svyvar}}
	\code{\link{lavaan}}
}
\examples{
###### A single group example #######

data(ess.dk)

dk.model <- "
  socialTrust ~ 1 + systemTrust + fearCrime
  systemTrust ~ 1 + socialTrust + efficacy
  socialTrust ~~ systemTrust
"
lavaan.fit <- lavaan(dk.model, data=ess.dk, auto.var=TRUE, estimator="MLM")
summary(lavaan.fit)
survey.design <- svydesign(ids=~intnum, data=ess.dk)

survey.fit <- lavaan.survey(lavaan.fit=lavaan.fit, survey.design=survey.design)
summary(survey.fit)

attr(survey.fit,"creff.svy")

# A test for R CMD CHECK 
stopifnot(abs(attr(survey.fit,"creff.svy")[2] - 1.2713339) < 1e-6)	
	
	
###### A multiple group example #######

data(HolzingerSwineford1939)

# The Holzinger and Swineford (1939) example - some model with complex restrictions
HS.model <- ' visual  =~ x1 + x2 + c(lam31, lam31)*x3
              textual =~ x4 + x5 + c(lam62, lam62)*x6
              speed   =~ x7 + x8 + c(lam93, lam93)*x9 
             speed ~ textual 
             textual ~ visual'

# Fit multiple group per school
fit <- lavaan(HS.model, data=HolzingerSwineford1939,
              auto.var=TRUE, auto.fix.first=TRUE, group="school",
              auto.cov.lv.x=TRUE, estimator="MLM")
summary(fit, fit.measures=TRUE)

# Create fictional clusters in the HS data
set.seed(20121025)
HolzingerSwineford1939$clus <- sample(1:100, size=nrow(HolzingerSwineford1939), replace=TRUE)
survey.design <- svydesign(ids=~clus, data=HolzingerSwineford1939)

summary(fit.survey <- lavaan.survey(fit, survey.design))

# Obtain a "relative efficiency" measure:
attr(fit.survey, "creff.svy")

#TODO: stopifnot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survey}
\keyword{models}
\keyword{regression}
\keyword{robust}
\keyword{multivariate}
