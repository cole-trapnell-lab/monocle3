#include <Rcpp.h>
using namespace Rcpp;

// Compute jaccard coefficient between nearest-neighbor sets
//
// Weights of both i->j and j->i are recorded if they have intersection. In this case
// w(i->j) should be equal to w(j->i). In some case i->j has weights while j<-i has no
// intersections, only w(i->j) is recorded. This is determinded in code `if(u>0)`.
// The original method described in the phenograph paper is used to calculate the weight.
//
// Author: Chen Hao, Date: 25/09/2015; updated by Xiaojie Qiu Nov. 12, 2017

NumericMatrix jaccard_coeff_cpp(NumericMatrix idx, bool weight) {
  int nrow = idx.nrow(), ncol = idx.ncol(), r = 0;
  NumericMatrix weights(nrow*ncol, 3);

  for(int i = 0; i < nrow; i ++) {
    for(int j = 0; j < ncol; j ++) {
      int k = idx(i,j) - 1;

      weights(r, 0) = i + 1;
      weights(r, 1) = k + 1;
      weights(r, 2) = 1;

      if(weight == TRUE) {

        NumericVector nodei = idx(i, _);
        NumericVector nodej = idx(k, _);

        int u = intersect(nodei, nodej).size();  // count intersection number
        int v = 2 * ncol - u;  // count union number

        if(u>0) {
          // weights(r, 0) = i + 1;
          // weights(r, 1) = k + 1;
          // weights(r, 2) = u / (2.0 * ncol - u) / 2;  // symmetrize the graph

          weights(r, 2) = (double) u / (double) v;  // normalize the values
        }
      }

      r ++;

    }
  }

  weights(_, 2) = weights(_, 2) / max(weights(_, 2));

  return weights;
}

// [[Rcpp::export]]
NumericMatrix jaccard_coeff(SEXP R_idx, SEXP R_weight) {
  NumericMatrix idx(R_idx);
  bool weight = as<bool>(R_weight);

  return jaccard_coeff_cpp(idx, weight);
}

NumericMatrix pnorm_over_mat_cpp(NumericMatrix num_links_ij, NumericMatrix var_null_num_links) {
  int n = num_links_ij.nrow();
  NumericMatrix tmp(n, n);

  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j ++) {
      // tmp(i, j) = Rcpp::pnorm( num_links_ij(i, j), 0.0, sqrt(var_null_num_links(i, j)), bool lower = false, bool log = false );
      tmp(i, j) = R::pnorm(num_links_ij(i, j), 0.0, sqrt(var_null_num_links(i, j)), 0, 0);
    }
  }
  return tmp;
}

// [[Rcpp::export]]
NumericMatrix pnorm_over_mat(SEXP R_num_links_ij, SEXP R_var_null_num_links) {
  NumericMatrix num_links_ij(R_num_links_ij);
  NumericMatrix var_null_num_links(R_var_null_num_links);

  return pnorm_over_mat_cpp(num_links_ij, var_null_num_links);
}


/***
edges$C = jaccard_dist
edges = subset(edges, C != 0)
edges$C = edges$C/max(edges$C)

*/
