# self training
self.train.kernel <- function (K, y, type = 'response', C = 1, 
                               cache = 40, tol = 0.001,
                               shrinking = TRUE, thrConf = 0.9,
                               maxIts = 10, percFull = 1, verbose = FALSE) 
{
  if(dim(K)[1]!=dim(K)[2])
    stop(" kernel matrix not square!")
  
  learner = "ksvm"
  learner.pars = list(type = "C-svc", prob.model = TRUE,
                      kernel = "matrix", C = C, cache = cache,
                      tol = tol, shrinking = shrinking)
  pred <- function(m,d)
  {
    p <- predict(m, as.kernelMatrix(d), type="prob")
    data.frame(cl = colnames(p)[apply(p, 1, which.max)],
               p = apply(p, 1, max))
  }
  pred.pars = list()
  
  data <- data.frame(y, K)
  N <- nrow(data)
  it <- 0
  sup <- which(!is.na(y))
  repeat
  {
    it <- it + 1
    model <- do.call(learner, c(list(K[sup, sup], data[sup, 1]), 
                                learner.pars))
    test.kernel <- as.kernelMatrix(K[,sup][, SVindex(model), drop = FALSE])
    pred.type <- predict(model, test.kernel, type = type)
    probPreds <- do.call(pred, c(list(model, test.kernel), pred.pars))
    new <- which(probPreds[, 2] > thrConf)
    if (verbose) 
      cat("IT.", it, "\t nr. added exs. =", length(new), "\n")
    if (length(new))
    {
      index <- setdiff(new, sup)
      sup <- unique(c(sup, new))
      data[index, 1] <- probPreds[index, 1]
    }
    else break
    if (it == maxIts || length(sup)/N >= percFull) 
      break
  }
  return(pred.type)
}