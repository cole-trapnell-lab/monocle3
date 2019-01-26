#' Function to automatically learn the structure of data by either using L1-graph or the spanning-tree formulization
#' @param X the input data DxN
#' @param C0 the initialization of centroids
#' @param G graph matrix with side information where cannot-link pair is 0
#' @param maxiter maximum number of iteraction
#' @param eps relative objective difference
#' @param gstruct graph structure to learn, either L1-graph or the span-tree
#' @param L1.lambda regularization parameter for inverse graph embedding
#' @param L1.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @param L1.sigma bandwidth parameter
#' @param nn number of nearest neighbors
#' @param L1.timeout a positive integer value specifying the number of seconds after which a timeout will occur. If zero, then no timeout will occur. (This is a parameter passed to lp.control function)
#' @param verbose emit results from iteraction
#' @return a list of X, C, W, P, objs
#' X is the input data
#' C is the centers for principal graph
#' W is the pricipal graph matrix
#' P is the cluster assignment matrix
#' objs is the objective value for the function
calc_principal_graph <- function(X, C0, G,
                            maxiter = 10, eps = 1e-5,
                            gstruct = c('l1-graph', 'span-tree'),
                            L1.lambda = 1,
                            L1.gamma = 0.5,
                            L1.sigma = 0.01,
                            nn = 5,
                            L1.timeout = 1800,
                            verbose = T) {

  C <- C0;
  D <- nrow(C); K <- ncol(C)

  # construct only once for all iterations
  if(gstruct =='l1-graph'){

    # low triangular sum_i sum_{j < i}
    G_tmp <- G
    G_tmp[upper.tri(G_tmp)] <- 0
    row_col <- which(G_tmp == 1, arr.ind = T) #lower triangle
    row <- row_col[, 1]; col <- row_col[, 2]
    nw <- length(row) #number of non-zero elements
    nvar <- nw + K*D #total number of variable; K*D: dimension of C

    rc <- new.env(hash=T, parent=emptyenv())
    # rc = java.util.HashMap; %hash-map: change the low triange into a vector and find the matrix index# change to a vector for linear programing
    for(i in 1:nw) {
      key_ij <- as.character(row[i] + col[i]*K) #index for the W_ij for the vector form of the weight matrix
      rc[[key_ij]] <- i
    }

    # construct A and b
    A <- matrix(nrow = 2 * D * K, ncol = nvar) #sparseMatrix(2 * D * K, nvar) # (i, j, x = x))
    b <- matrix(nrow = 2 * D * K, ncol = 1) #sparseMatrix(2 * D * K, 1) #dimension is (2D * K) X 1

    for(i in 1:K){
      # print(i)
      nn_i <- which(G[,i]==1)
      a <- matrix(0, nrow = 2 * D, ncol = nvar)
      if (length(nn_i) >= 1){
        for(jj in 1:length(nn_i)){
          j <- nn_i[jj]
          key_ij <- as.character(i+j*K)
          if(i < j) {
            key_ij <- as.character(j + i*K) #index for the neighbor for the vector form of the graph
          }
          pos_ij <- rc[[key_ij]]
          #put this into the corresponding column
          a[, pos_ij] <- c(-X[, j], X[, j]) # -a < x < a; eq. 17. i-th sample's constraints. each sample has D genes
        }
      }

      start_i <- nw + (i-1)*D + 1
      end_i <- start_i + D-1 #why do we need D - 1 here? should be just D, right?
      #those are the columns corresponding to the episilon
      a[, start_i:end_i] = -rbind(matlab::eye(D,D), matlab::eye(D,D)) # |a| < b =>  a - b <= 0 & -b -a >= 0; eye: comes from the b (I matrix)
      A[((i - 1) * 2 * D + 1):((i) * 2 * D), ] <- a
      b[((i - 1) * 2 * D + 1):((i) * 2 * D), ] <- c(-X[, i], X[, i])
    }
  }

  objs <- c()
  lp_vars <- c()
  for(iter in 1:maxiter){

    norm_sq <- repmat(t(colSums(C^2)), K, 1) #this part calculates the cost matrix Phi
    Phi <- norm_sq + t(norm_sq) - 2 * t(C) %*% C
    if(gstruct == 'l1-graph'){
      val <- matrix(0, nrow = nw, ncol = 1)
      for(i in 1:nw) {
        val[i] = Phi[row[i], col[i]]
      }

      f <- matrix(c(2*val, L1.lambda * rep(1, K*D)), ncol = 1)

      # MATLAB solver
      #options = optimset( 'Display', 'off','Algorithm','interior-point');
      # [w_eta, obj_W] = linprog(f, A, b, [], [], [zeros(nw, 1); -Inf*matrix(K*D,nrow = 1)], [], lp_vars, options);

      # lpSolve package
      # prod.sol <- lp("min", f, A, constr.dir, b, compute.sens = TRUE)
      # obj_W <- prod.sol$obj.val #objective function value
      # w_eta <- prod.sol$solution

      #another approach:
      # nrow a nonnegative integer value specifying the number of constaints in the linear program.
      # ncol a nonnegative integer value specifying the number of decision variables in the linear program.
      lprec <- lpSolveAPI::make.lp(length(b), length(f), verbose="important")
      lpSolveAPI::set.objfn(lprec, f)
      for(i in 1:nrow(A)) {
        lpSolveAPI::add.constraint(lprec, A[i, ], "<=", b[i, ])
      }
      # for(j in 1:nrow()) {
      # 	set.bounds(lprec, lower = c(rep(0, nw), -Inf*rep(1, K*D)), columns = 1:length(b))
      # }
      lpSolveAPI::set.bounds(lprec, lower = c(rep(0, nw), -Inf*rep(1, K*D))) #				set.bounds(lprec, lower = c(rep(0, nw), -Inf*rep(1, K*D)), columns = 1:length(b))


      #set objective direction
      #lp.control(lpmodel,sense='max')

      #solve the model, if this return 0 an optimal solution is found
      lpSolveAPI::lp.control(lprec, timeout = L1.timeout, presolve = "rows")
      solve(lprec)

      obj_W <- lpSolveAPI::get.objective(lprec)
      w_eta <- lpSolveAPI::get.variables(lprec)

      # %                                                           lb(w, )                      ub
      # %             % Mosek solver
      # %             prob.c = f; prob.a = A;
      # %             prob.buc = b;
      # %             prob.blx = sparse( [zeros(nw, 1); -Inf.*ones(K*D,1)] );
      # %             [r,res] = mosekopt('minimize echo(0)',prob);
      # %             w_eta = res.sol.bas.xx;
      # %             obj_W = f'*w_eta;

      # recover results
      w <- w_eta[1:nw]
      W_tril <- as.matrix(sparseMatrix(row, col, x = w, dims = c(K, K))) #W_tril(row(i), row(i)) = w(i) and K should be dimension of the matrix
      W <- W_tril + t(W_tril)

      # warm start #it seems that the matlab version ignored warm start with the interior-point method
      lp_vars <- w_eta
    }
    else if(gstruct == 'span-tree'){
      ##########################use mst from igraph: ##########################
      g <- igraph::graph.adjacency(Phi, mode = 'lower', diag = T, weighted = T)
      g_mst <- igraph::mst(g)
      stree <- igraph::get.adjacency(g_mst, attr = 'weight', type = 'lower')
      stree_ori <- stree

      #convert to matrix:
      stree <- as.matrix(stree)
      stree <- stree + t(stree)

      W <- stree != 0
      obj_W <- sum(sum(stree))
    }
    else
      warning('graph structure %s is not supported yet.',gstruct)

    res = soft_assignment(X, C, L1.sigma)
    P <- res$P
    obj_P <- res$obj

    obj <- obj_W + L1.gamma * obj_P
    objs = c(objs, obj)
    if(verbose)
      message('iter = ', iter, ' obj = ', obj)

    if(iter > 1){
      relative_diff = abs( objs[iter-1] - obj) / abs(objs[iter-1]);
      if(relative_diff < eps){
        if(verbose)
          message('eps = ', relative_diff, ', converge.')
        break
      }
      if(iter >= maxiter){
        if(verbose)
          message('eps = ', relative_diff, ' reach maxiter.')
      }
    }

    C <- generate_centers(X, W, P, L1.gamma)

  }

  return(list(X = X, C = C, W = W, P = P, objs = objs))
}


#' function to reproduce the behavior of repmat function in matlab to replicate and tile an matrix
#' @param X matrix for tiling and replicate the data
#' @param m a numeric value for tiling a matrix
#' @param n a numeric value for tiling a matrix
#' @return a matrix
repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}



#' Function to calculate the third term in the objective function
#' @param X input data
#' @param C center of grap (D * K)
#' @param sigma bandwidth parameter
#' @return a matrix with diagonal element as 1 while other elements as zero (eye matrix)
soft_assignment <- function(X, C, sigma){

  D <- nrow(X); N <- ncol(X)
  K <- ncol(C)
  norm_X_sq <- repmat(t(t(colSums(X^2))), 1, K);
  norm_C_sq <- repmat(t(colSums(C^2)), N, 1);
  dist_XC <- norm_X_sq + norm_C_sq - 2 * t(X) %*% C

  # %% handle numerical problems 0/0 for P
  min_dist <- apply(dist_XC, 1, min) #rowMin(dist_XC)

  dist_XC <- dist_XC - repmat(t(t(min_dist)), 1, K )
  Phi_XC <- exp(- dist_XC / sigma)
  P <- Phi_XC / repmat(t(t(rowSums(Phi_XC))), 1, K)

  obj <- - sigma * sum( log( rowSums( exp(- dist_XC/sigma)) ) #why not \sigma * log (p_{ij})
                        - min_dist/ sigma );

  return(list(P = P, obj = obj))
}

#' Function to reproduce the behavior of eye function in matlab
#' @param X input data
#' @param W the pricipal graph matrix
#' @param P the cluster assignment matrix
#' @param param.gamma regularization parameter for k-means (the prefix of 'param' is used to avoid name collision with gamma)
#' @return A matrix C for the centers for principal graph
#'
generate_centers <- function(X, W, P, param.gamma){
  D <- nrow(X); N <- nrow(X)
  K <- ncol(W)
  # prevent singular
  Q <- 2 *( diag(colSums(W)) - W ) + param.gamma * diag(colSums(P))
  B <-  param.gamma * X %*% P;
  C <- B %*% solve(Q)   #equation 22

  return(C)
}



