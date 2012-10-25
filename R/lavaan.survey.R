
lavaan.survey <-
function(lavaan.fit, survey.design) {
  ov.names <- lavaan.fit@Data@ov.names[[1]]
  Dplus <- ginv(lavaan::duplicationMatrix(length(ov.names)))
  ov.formula <- as.formula(paste("~",paste(ov.names, collapse="+")))
  sample.nobs <- unlist(lavaan.fit@Data@nobs)
  
  Gamma <- vector("list", lavaan.fit@Data@ngroups)
  sample.cov <- vector("list", lavaan.fit@Data@ngroups)
  sample.mean <- vector("list", lavaan.fit@Data@ngroups)
  
  for(g in seq(lavaan.fit@Data@ngroups)) {
    if(lavaan.fit@Data@ngroups > 1) {
      survey.design.g <- subset(survey.design, eval(parse(text=sprintf("%s == '%s'", 
                              lavaan.fit@call$group, lavaan.fit@Data@group.label[[g]]))))
    } else { survey.design.g <- survey.design  }
    
    sample.cov.g <- as.matrix(svyvar(ov.formula, design=survey.design.g, na.rm=TRUE))  
    Gamma.cov.g <- attr(sample.cov.g, "var")
    Gamma.cov.g <- Dplus %*% Gamma.cov.g %*% t(Dplus)
    attr(sample.cov.g, "var") <- NULL
  
    # Check positive-definiteness of Gamma for the covariances
    stopifnot(sum(eigen(Gamma.cov.g)$values>1e-6) == ncol(Gamma.cov.g))
  
    sample.mean.g <- svymean(ov.formula, design=survey.design.g, na.rm=TRUE)  
    Gamma.mean.g <- attr(sample.mean.g, "var")
    Gamma.g <- as.matrix(Matrix::bdiag(Gamma.mean.g, Gamma.cov.g)) # TODO add offdiag
    attr(sample.mean.g, "var") <- NULL
    # Check positive-definiteness of Gamma for the means
    stopifnot(sum(eigen(Gamma.mean.g)$values>1e-6) == ncol(Gamma.mean.g))
  
    Gamma.g <- Gamma.g * sample.nobs[g]
    
    Gamma[[g]] <- Gamma.g
    sample.cov[[g]] <- sample.cov.g
    sample.mean[[g]] <- sample.mean.g
  }
  
  new.call <- lavaan.fit@call
  new.call$data <- NULL                # Remove any data argument
  new.call$sample.cov <- sample.cov    # Set survey covariances
  new.call$sample.mean <- sample.mean  # Set survey means
  new.call$sample.nobs <- sample.nobs  
  new.call$estimator <- "MLM"  # Always use Satorra-Bentler method
  new.call$NACOV <- Gamma      # Set asy covariance matrix of sample means and covariances
  
  new.fit <- eval(new.call) # Run lavaan with the new arguments
  
  attr(new.fit, "creff.svydesign") <-
    tryCatch(sqrt(diag(vcov(new.fit)) / diag(vcov(lavaan.fit))),
    error = function(e) {warning("Attempt to compute vcov failed for one of the models."); NULL})
  
  new.fit
}
