# inferring functionally related proteins using networks
net.infer <- function(target, kernel, top=NULL, cross=0,
                      C = 1, nu = 0.2, epsilon = 0.1, cache1 = 40,
                      tol1 = 0.001, shrinking1 = TRUE, cache2 = 40,
                      tol2 = 0.001, shrinking2 = TRUE)
{
  # find new list
  node <- rownames(kernel)
  n <- length(node)
  if(length(target)<=1) stop("size of list is too small")
  index <- match(target,node)!='NA'
  new.target <- target[!is.na(index)]
  n1 <- length(na.omit(index))
  if(n-2*n1<=0) stop("size of list is too large")
  
  
  #### OCSVM
  
  # train
  new.index <- !is.na(match(node,new.target))
  train.kernel <- kernel[new.index,new.index]
  model <- ksvm(train.kernel, type="one-svc",kernel="matrix",
                nu = nu, epsilon = epsilon, cache = cache1,
                tol = tol1, shrinking = shrinking1)
  # predict
  test.kernel <- as.kernelMatrix(kernel[,new.index][,SVindex(model), drop=FALSE]) 
  pred.decision <- predict(model, test.kernel,type='decision')
  
  
  # assign class 1
  y <- ifelse(new.index==TRUE,1,'NA')
  
  pred <- cbind(node,pred.decision,y)
  sort.pred <- pred[order(as.numeric(pred[,2])),]
  sort.pred2 <- subset(sort.pred,sort.pred[,3]=='NA')
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[,2])),]
  target0 <- as.vector(new.sort.pred2[1:n1,1])
  
  # assign pseudo class 0
  new.index <- !is.na(match(node,target0))
  y <- ifelse(new.index==TRUE,0,y)
  
  
  
  #### SVM
  
  # train
  model <- ksvm(kernel[y!="NA", y!="NA"], as.numeric(y[y!="NA"]),
                type = "C-svc", C = C, kernel="matrix", cross=cross,
                cache = cache2, tol = tol2, shrinking = shrinking2)
  
  
  # predict
  test.kernel <- as.kernelMatrix(kernel[, y!="NA"][,SVindex(model), drop=FALSE]) 
  pred.decision <- predict(model, test.kernel,type='decision')
  
  
  
  
  pred <- cbind(node,pred.decision,y)
  colnames(pred)[2] <- 'decision values'
  
  sort.pred <- pred[order(as.numeric(pred[,2]),decreasing = TRUE),]
  sort.pred2 <- subset(sort.pred,sort.pred[,3]=='NA')
  new.sort.pred2 <- sort.pred2[order(as.numeric(sort.pred2[,2]),
                                     decreasing = TRUE),]
  
  
  if( is.null(top) ) {top=n-2*n1}  
  if( (n-2*n1)<top ) {top=n-2*n1}
  pred.svm <- list()
  pred.svm$list <- new.target
  pred.svm$error <- error(model)
  if (cross > 0) {pred.svm$CVerror=cross(model)}
  pred.svm$top <- as.vector(new.sort.pred2[1:top,1])
  pred.svm$score <- as.numeric(new.sort.pred2[1:top,2])
  
  pred.svm
}