\name{liss}
\alias{liss}
\docType{data}
\title{
  Internet use in a panel data set.
}
\description{
The longitudinal internet studies for the social sciences (LISS) panel is a web survey panel recruited by probability sampling of households. All household members participate in the survey. This dataset contains four waves of data (2008-2011) on panel participants' internet usage.
}
\usage{data(liss)}
\format{
  A data frame with 7369 observations on the following 6 variables.
  
  \describe{
    \item{\code{nohouse_encr}}{Household identifier.}
    \item{\code{nomem_encr}}{Person identifier.}
    \item{\code{cs08a247}}{Natural logarithm of number of hours internet usage at home per week in 2008.}
    \item{\code{cs09b247}}{Natural logarithm of number of hours internet usage at home per week in 2009.}
    \item{\code{cs10c247}}{Natural logarithm of number of hours internet usage at home per week in 2010.}
    \item{\code{cs11d247}}{Natural logarithm of number of hours internet usage at home per week in 2011.}
  }
}
\source{
  Data were obtained from http://www.lissdata.nl/dataarchive/. See also http://www.lissdata.nl/dataarchive/data_variables/view/795. 
}
\references{
Scherpenzeel A.C. (2011). "Data Collection in a Probability-Based Internet Panel: 
  How the LISS Panel was Built and How it Can be Used." 
  Bulletin of Sociological Methodology, 109(1), 56-61.

Oberski, D.L. (2014). lavaan.survey: An R Package for Complex Survey Analysis
  of Structural Equation Models. Journal of Statistical Software, 57(1), 1-27.
  \url{http://www.jstatsoft.org/v57/i01/}.
}
\examples{
data(liss)

# Estimating the reliability of internet usage with the "quasi-simplex" 
#   (a.k.a. Gaussian latent Markov) model. 

# A quasi-simplex model for four time points, setting the error variance 
#   to be estimated equal over time.
model.liss <- "
   cs08 =~ 1 * cs08a247
   cs09 =~ 1 * cs09b247
   cs10 =~ 1 * cs10c247
   cs11 =~ 1 * cs11d247
 
   cs09 ~ cs08
   cs10 ~ cs09
   cs11 ~ cs10
 
   cs08a247 ~~ vare * cs08a247
   cs09b247 ~~ vare * cs09b247
   cs10c247 ~~ vare * cs10c247
   cs11d247 ~~ vare * cs11d247
 
   cs08 ~~ vart08 * cs08
 
   reliab.ratio := vart08 / (vart08 + vare)
"

# Fit the model using listwise deletion
fit.liss <- lavaan(model.liss, auto.var = TRUE, meanstructure = TRUE, 
  int.ov.free = TRUE, data = liss)
  
# Fit the model accounting for nesting of respondents within households
des.liss <- svydesign(ids = ~nohouse_encr, prob = ~1, data = liss)
fit.liss.surv <- lavaan.survey(fit.liss, des.liss)
fit.liss.surv

# Complex survey inference on the reliability of interest:
parameterEstimates(fit.liss.surv)[24, ]

## To deal with missing data (including attrition), multiple imputation can be used.
## For example using the mice library (although any MI software is suitable)

## Uncomment below to run this time-intensive analysis
## NOT RUN:

# set.seed(20140221)
# library("mice") 
# liss.imp <- mice(liss, m = 100, method = "norm", maxit = 100)

## Turn the mice object into a list() of imputed datasets
# liss.implist <- lapply(seq(liss.imp$m), function(im) complete(liss.imp, im))

## After obtaining the list of imputed datasets, 
##  use the mitools package to turn it into an imputation list
# library("mitools")
# liss.implist <- imputationList(liss.implist)

## Give the imputation list as data to a svydesign object
# des.liss.imp <- svydesign(ids = ~nohouse_encr, prob = ~1, data = liss.implist)

## lavaan.survey can be used as usual, using the 
##    svydesign object that has an imputation list as data
## Standard errors and chi-square tests will account for both the clustering and the 
##   imputation uncertainty applying Rubin's rules. 
# fit.liss.surv.mi <- lavaan.survey(fit.liss, des.liss.imp)
# fit.liss.surv.mi

## After this analysis, we can again perform inference on the reliability of interest:
# parameterEstimates(fit.liss.surv.mi)[24, ]

}
\author{
  Daniel Oberski - \url{http://daob.nl/} - \email{daniel.oberski@gmail.com}
}
\keyword{datasets}
