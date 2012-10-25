subset.survey.design <-
function(survey.design, subset.call) {
    log <- eval(substitute(subset.call), envir=survey.design$variables, 
                 enclos = parent.frame())
    survey.design$cluster <- subset(survey.design$cluster, log)
    survey.design$strata <- subset(survey.design$strata, log)
    survey.design$prob <- survey.design$prob[log]
    survey.design$allprob <- subset(survey.design$allprob, log)
    survey.design$variables <- subset(survey.design$variables, log)
    survey.design$fpc$sampsize <- subset(survey.design$fpc$sampsize, log)
    
    survey.design
}
