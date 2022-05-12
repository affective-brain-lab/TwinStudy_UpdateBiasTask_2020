# control for confounds
# indep is the independent variables
# x is datamatrix, ctrl is names of variables to be controlled
# ctrl.factor is boolean array shpwing which of the variables to be controlled are factors
control_confounds <- function(x, indep, ctrl, ctrl.factor=rep(FALSE, length(ctrl))) {
  
  stopifnot(length(ctrl)==length(ctrl.factor))
  
  for(i in which(ctrl.factor)) {
    x[[ctrl[i]]] <- as.factor(x[[ctrl[i]]])
  }
  
  ctrl <- paste(ctrl,collapse="+")
  cat("controlling for:",ctrl,"\n")    
  txt <- paste(" res <- lm(",indep,"~",ctrl,",data=x, na.action=na.exclude)")
  eval(parse(text=txt))
  resid(res)
}
