#' Robust Convex Clustering
#' 
#' \code{robustcvxclust_path} estimates the robust convex clustering path using ADMM.
#' Required inputs include a data matrix \code{X} (rows are the observations; columns are the features), a vector of weights \code{wt}, and a sequence of regularization parameters \code{lambda}. The penalty norm supported is 2-norm.
#' 
#' @param X The data matrix to be clustered. The rows are the observations, and the columns are the features.
#' @param wt A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids.
#' @param lambda A sequence of regularization parameters.
#' @param rho Augmented Lagrangian penalty parameter.
#' @param tau The robustification parameter in huber loss.
#' @param tol_abs The convergence tolerance (absolute).
#' @param max_iter The maximum number of iterations.
#'
#' @return \code{U} A list of centroid matrices.
#' @return \code{W} A list of centroid matrices.
#' @return \code{V} A list of centroid difference matrices.
#' @return \code{Y} A list of Lagrange multiplier matrices.
#' @return \code{Z} A list of Lagrange multiplier matrices.
#' @export
#' @author Chenyu Liu, Qiang Sun, Kean Ming Tan
#' @useDynLib Rcvxclustr
#' @examples 
#' ## Create random problems
#' p <- 10
#' n <- 20
#' seed <- 1234
#' mu1 <- rnorm(p,0,1)
#' X1 <- mvrnorm((n/2),mu=mu1,Sigma=diag(1,p))
#' mu2 <- c(rnorm(p/2,3,1),rnorm(p/2,-3,1))
#' X2 <- mvrnorm((n/2),mu=mu2,Sigma=diag(1,p))
#' X <- rbind(X1,X2)
#' out.percent <- 0.06
#' n.out <- as.integer(n*out.percent)
#' nth <- sort(sample.int(n.out,1:n))
#' for(i in nth){
#'  outliers <- sample(runif(n=p*0.2,min = 10,max = 20))
#'  X[i,sample(p)[1:(p*0.2)]] <- outliers
#' }
#' wt <- rep(1,n*(n+1)/2)
#' rho <- 1
#' tau <- 3
#' lambda <- seq(1,5,1)
#' sol <- robustcvxclust_path(X,wt,i,rho,tau,tol_abs,max_iter)
robustcvxclust_path <-
  function(X,wt,lambda,rho,tau,tol_abs=1e-05, max_iter=1e5){
    call <- match.call()
    nLambda <- length(lambda)
    n <- nrow(X)
    p <- ncol(X)
    E <- create_E_matrix(n)$E
    nK <- nrow(E)
    list_U <- vector(mode="list",length=nLambda)
    list_W <- vector(mode="list",length=nLambda)
    list_V <- vector(mode="list",length=nLambda)
    list_Y <- vector(mode="list",length=nLambda)
    list_Z <- vector(mode="list",length=nLambda)
    iter_vec <- integer(nLambda)
    for(id in 1:nLambda){
      lam <- lambda[id]
      rcc <- robustcvxclust(X,W,V,Y,Z,max_iter=max_iter,rho=rho,tau=tau,lambda=lam,wt=wt,tol_abs=tol_abs)
      iter_vec[id] <- rcc$iter
      list_U[[id]] <- rcc$U
      list_W[[id]] <- rcc$W
      list_V[[id]] <- rcc$V
      list_Y[[id]] <- rcc$Y
      list_Z[[id]] <- rcc$Z
    }
    robustcvxclust_obj <- list(U=list_U,W=list_W,V=list_V,Y=list_Y,Z=list_Z,nLambda=id,
                               iters=iter_vec,call=call)
    class(robustcvxclust_obj) <- "robustcvxclustobject"
    return(robustcvxclust_obj)
  }
