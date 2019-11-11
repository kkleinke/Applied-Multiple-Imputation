test.pan <- function(y, x, xcol, zcol, group, paniter = 500, intercept=TRUE, m=5){
  # arguments:
  ############
  # y: matrix of incomplete response variables (variables to be imputed)
  # x: matrix of complete predictors
  # xcol: vector of integers indicating which columns of x are fixed effects
  # zcol: vector of integers indicating which columns of x are random effects
  # group: group identifier
  # paniter: number of iterations before each imputation is drawn
  # intercept: shall intercept be included in model?
  # m: number of multiple imputations to be created
  #
  #output: list of length m of completed data sets
  ################################################################
  
  require(pan)
  if (intercept) {
    x <- cbind(1, data.matrix(x))
  }
  ry=!is.na(y)
  subj <- match(group, unique(group))
  sortgroups <- any(diff(subj) < 0)
  if (sortgroups) {
    dfr <- data.frame(group = group, ry = ry, index = seq(1,length(ry)))
    dfr <- dfr[order(dfr$group), ]
    group <- group[dfr$index]
    y <- y[dfr$index]
    x <- x[dfr$index, ]
    ry <- ry[dfr$index]
    subj <- subj[dfr$index]
  }
  y1 <- matrix(as.numeric(y), ncol = 1)
  y1[!ry, 1] <- NA
  
  prior <- list(a = ncol(y1), Binv = diag(rep(1, ncol(y1))),
                c = ncol(y1) * length(zcol), Dinv = diag(rep(1, ncol(y1) *
                                                               length(zcol))))
  if (length(subj) != nrow(y1))
    stop("No class variable")
  
  s1 <- round(runif(1, 1, 10^7))
  result <- pan(y1, subj, x, xcol, zcol, prior, seed = s1,
                iter = paniter)
  return(result)
}