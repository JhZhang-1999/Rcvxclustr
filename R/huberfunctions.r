##### Huber Loss #####

### Main Function ###
#' @export
Huber_ADMM <- function(X, initV=NULL,initY=NULL,initW=NULL,initZ=NULL,
                       maxit=1000,rho=1,tau=3,lam=0.5,eps=1e-3,wt=NULL){
  n = nrow(X)
  p = ncol(X)
  I = diag(rep(1,n))
  ### contribute E
  for(i in seq(from=n-1,to=1,by=-1)){
    if(i>1){
      temp=diag(rep(-1,i))
    }else{temp=-1}
    temp1=cbind(rep(1,i),temp)
    if(n-i-1==0){
      E=temp1
    }else{E=rbind(E,cbind(matrix(0,ncol=n-i-1, nrow=i),temp1))}
  }
  k = dim(E)[1] # the amount of E's rows

  ### Generate initial V,Y,W,Z, if no initial values then generate randomly
  if(is.null(initV)==TRUE){
    V <- matrix(rnorm(k*p),nrow = k,ncol = p)
  }else{V <- initV}
  if(is.null(initY)==TRUE){
    Y <- matrix(rnorm(k*p),nrow = k,ncol = p)
  }else{Y <- initY}
  if(is.null(initW)==TRUE){
    W <- matrix(rnorm(n*p),nrow = n,ncol = p)
  }else{W <- initW}
  if(is.null(initZ)==TRUE){
    Z <- matrix(rnorm(n*p),nrow = n,ncol = p)
  }else{Z <- initZ}
  if(is.null(wt)==TRUE){
    wt=rep(1,k)
  }

  sigma <- lam*wt/rho

  A <- matrix(tau,nrow=n,ncol=p) # A is tau formed as a matrix

  for(t in 2:maxit){
    Vold <- V
    Yold <- Y
    Wold <- W
    Zold <- Z
    ### Update U
    U <- solve(t(E)%*%E+I, t(E)%*%(Vold+Yold)+Wold+Zold) # A1.3.(a)
    ### Update W
    B <- abs(rho*(X-U+Zold))/(1+rho)
    index <- which(B<=A,arr.ind=T)
    W <- X - Soft(X-U+Zold,tau/rho) # A1.3.(b)
    W[index] <- (X[index]+rho*(U[index]-Zold[index]))/(1+rho) # A1.3.(b)
    ### Update V
    v <- E%*%U - Yold # v is e in the paper
    v.norm <- apply(v,1,function(x){norm(x,type = "2")}) # the p2 norms of the rows
    V <- v*positive_part(1-sigma/v.norm) # A1.3.(c)
    ### Update Y
    Y <- Yold-rho*(E%*%U-V) # A1.3.(d)
    ### Update Z
    Z <- Zold-rho*(U-W) # A1.3.(e)
    if(norm(W-Wold,type="2")<eps){break} # if tolerance level met then stop iteration
    Uold <- U
  }
  return(list(W=W,V=V,iteration=t,convergence=norm(W-Wold,type="2")))
}
#####################################################

### Soft-thresholding function ###
#' @export
Soft <- function(a,b){
  if(b<0) stop("Can soft-threshold by a nonnegative quantity only.")
  return(sign(a)*pmax(0,abs(a)-b))
}

### Positive Part function ###
#' @export
positive_part <- function(x){
  a <- x
  for(i in 1:length(x)){
    if(x[i]<0){a[i]=0}
  }
  return(a)
}

### The Proposed Weight Function ###
#' @export
weight_delta <- function(X, delta=15, zeta=0.1){
  n <- dim(X)[1]
  distance <- matrix(0,nrow=n,ncol=n)
  for(i in seq(from=1,to=n-1,by=1)){
    for(j in seq(from=i+1,to=n, by=1)){
      distance[i,j]<-norm(X[i,] - X[j,],type = "2")
      distance[j,i]<-norm(X[i,] - X[j,],type = "2")
    }
  }
  wt <- matrix(0,n,n) # n*n matrix, all elements=0
  for(i in 1:n){
    for(j in 1:n){
      if(distance[i,j]==0){
        wt[i,j]=0
      }
      else if(distance[i,j]<delta){
        wt[i,j]=exp(-0.5*zeta*distance[i,j])
      }
      else{wt[i,j]=exp(-0.5*zeta*delta)}
    }
  }
  return (wt)
}

