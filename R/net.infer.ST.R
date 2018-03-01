# inferring functionally related proteins with self training
net.infer.ST <- function(target, kernel, top = NULL, C = 1, nu = 0.2, 
                         epsilon = 0.1, cache1 = 40, tol1 = 0.001, shrinking1 = TRUE, 
                         cache2 = 40, tol2 = 0.001, shrinking2 = TRUE, thrConf = 0.9,
                         maxIts = 10, percFull = 1, verbose = FALSE) 
{
  node <- rownames(kernel)
  n <- length(node)
  if (length(target) <= 1) 
    stop("size of list is too small")
  index <- match(target, node) != "NA"
  new.target <- target[!is.na(index)]
  n1 <- length(na.omit(index))
  if (n - 2 * n1 <= 0) 
    stop("size of list is too large")
  
  
  #### OCSVM
  new.index <- !is.na(match(node, new.target))
  train.kernel <- kernel[new.index, new.index]
  model <- ksvm(train.kernel, type = "one-svc", kernel = "matrix", 
                C = C, nu = nu, epsilon = epsilon, cache = cache1, tol = tol1, 
                shrinking = shrinking1)
  test.kernel <- as.kernelMatrix(kernel[, new.index][, SVindex(model), 
                                                     drop = FALSE])
  pred.decision <- predict(model, test.kernel, type = "decision")
  y <- ifelse(new.index == TRUE, 1, "NA")
  pred <- cbind(node, pred.decision, y)
  sort.pred <- pred[order(as.numeric(pred[, 2])), ]
  sort.pred2 <- subset(sort.pred, sort.pred[, 3] == "NA")
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[, 
                                                           2])), ]
  target0 <- as.vector(new.sort.pred2[1:n1, 1])
  new.index <- !is.na(match(node, target0))
  y <- ifelse(new.index == TRUE, 0, y)
  
  # self training
  pred.decision  <- self.train.kernel(kernel, factor(y, levels=c(0,1)), type = 'decision',
                                      C = C, cache = cache2, tol = tol2, shrinking = shrinking2,
                                      thrConf = thrConf, maxIts = maxIts,
                                      percFull = percFull, verbose = verbose)
  
  pred <- cbind(node, pred.decision, y)
  colnames(pred)[2] <- "decision values"
  sort.pred <- pred[order(as.numeric(pred[, 2]), decreasing = TRUE), 
                    ]
  sort.pred2 <- subset(sort.pred, sort.pred[, 3] == "NA")
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[, 
                                                           2]), decreasing = TRUE), ]
  if (is.null(top)) {
    top = n - 2 * n1
  }
  if ((n - 2 * n1) < top) {
    top = n - 2 * n1
  }
  pred.svm <- list()
  pred.svm$list <- new.target
  pred.svm$error <- error(model)
  
  pred.svm$top <- as.vector(new.sort.pred2[1:top, 1])
  pred.svm$score <- as.numeric(new.sort.pred2[1:top, 2])
  pred.svm
}