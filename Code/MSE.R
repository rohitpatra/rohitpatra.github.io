#------------------------------------------------------------------------------------
#-----------Function to calculate the Maximum score estimator given a----------------
				# matrix of covariatesX and binary responses y-----------------------
#------------------------------------------------------------------------------------
# Author: Rohit Kumar Patra
# Date:  19 Jan 2016
# See: stat.ufl.edu/~rohitpatra/MSEcompute.html for more information about the code

library(Rcplex)
library(Matrix)

## Include beta.1 but make lower bound and upper bound 1
DefineControl <- function(X,y,w=NULL, search.width=10) {
	# See equation 3 of Florios and Skouras 2008
	if (!is.matrix(X)) X <- as.matrix(X)
	if (!is.vector(y)) y <- as.vector(y)
	n <- nrow(X)
	dd <- ncol(X)
	X.0 <- as.matrix(X[,1])
	X.1 <- as.matrix(X[,-1])
	if (is.null(w)) w <- rep(-1, n) # Weights for MSE is always 1.
	##-------------- The constant for objective function---------
	cvec<- c(w, rep(0, dd))

	#------------------------------------------------------------------------------------
	##-----The Matrix Amat and vector bvec for  for inequality constraint A* var \le b---
	#------------------------------------------------------------------------------------

	M <- abs(X.0)+ rowSums(abs(X.1))*search.width
	A.beta <- diag(1- 2*y)%*%X
	Amat <- cbind(diag(as.vector(M)), A.beta)
	Amat <- Matrix(Amat, sparse=TRUE)
	bvec <- as.vector(M)

	#------------------------------------------------------------------------------------
	##-----Upper and lower bound for the variables---------------------------------------
	#------------------------------------------------------------------------------------


	lb <- c(rep(0,n), 1,rep(-search.width, dd-1))
	# We set both the lower bound and upper for beta.1 to be 1.
	ub <- c(rep(1,n),1, rep(search.width, dd-1)) #

	##-----Variable types z1,..,zn are binary and beta.1,...., beta.dd are continuous.---------

	vtype <- c(rep("B", n), rep("C", dd))

	# The first n variable are z.1 through z.n and they are binary. The next dd  variables are beta. Beta.1 is fixed at one but the rest are in the window [-search.width, search.width]

 	### OUTPUT
	list(cvec=cvec,
		 Amat=Amat,
		 w=w,
		 lb=lb,
		 ub=ub,
		 bvec=bvec,
		 vtype=vtype,
		 dd=dd,
		 varlength=length(cvec) );
}
MSE <- function(X, y, search.width=10, tol= 0.02){
	#-----Evaluating the constraint variable for Cplex
	Cons <-  DefineControl(X, y, w=NULL, search.width=5)
	t <- Sys.time()
	out <- Rcplex(cvec=Cons$cvec, Amat=Cons$Amat, bvec=Cons$bvec,
					  Qmat=NULL, lb=Cons$lb, ub=Cons$ub,
					  objsense="min",sense="L", vtype=Cons$vtype,
					  control=list(tilim=3600,epgap=tol)) #,mipemphasis=3
	t1 <- Sys.time()
	#See  help of Rcplex for variable names.
	opt.beta <- out$xopt[-(1:nrow(X))] # The last d variable of optimizing vector.
	opt.beta <- opt.beta/norm(opt.beta, "2")
	list(out = out,
	 	 beta.hat = opt.beta,
	 	 Mscore = -out$obj,
	 	 time = t1-t);
}

#----------------------------------------------------------------------------
 							# Trial run codes
#----------------------------------------------------------------------------
# X <- cbind(c(2,-3,1),c(1,1,1))
# y <- c(0,1,0)
# ans <- MSE(X,y)
# Answer should be: ans$obj  = -2, ans$beta.hat = c(1, -10).

#----------------------------------------------------------------------------
							# Random data Application from the paper
#---------------------------------------------------------------------------
# n <- 200
# d <- 3
# X <- matrix(runif(n*d, -1,1), n,d)
# beta.0 <- rep(1,d)
# ind <- X%*%beta.0
# sdx <- (1 +  rowSums(X^2))^(-1)
# err <- rnorm(n,0,sdx)
# y <- as.vector((ind+err>0)*1)
# # plot(ind,y)
# ans <- MSE(X,y)

# print("#----------------------------------------------------------------------------")
# cat(" The MSE for this data set is:  " )
# print(ans$beta.hat, digits=3)
# print(paste("It took ", format(ans$time,  digits=2), " secs to compute it"))


