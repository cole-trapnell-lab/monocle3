#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
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

  // Temporary storage for sorting
  std::vector<double> v_i(ncol);
  std::vector<double> v_k(ncol);

  for(int i = 0; i < nrow; i ++) {
    // Fill and sort v_i once per row i
    if (weight) {
      for(int c = 0; c < ncol; c++) v_i[c] = idx(i, c);
      std::sort(v_i.begin(), v_i.end());
    }

    for(int j = 0; j < ncol; j ++) {
      int k = idx(i,j) - 1;

      weights(r, 0) = i + 1;
      weights(r, 1) = k + 1;
      weights(r, 2) = 1;

      if(weight == TRUE) {
        // Fill and sort v_k
        for(int c = 0; c < ncol; c++) v_k[c] = idx(k, c);
        std::sort(v_k.begin(), v_k.end());

        // Standard set_intersection logic to count overlaps
        int u = 0;
        int idx_i = 0;
        int idx_k = 0;
        
        while (idx_i < ncol && idx_k < ncol) {
            if (v_i[idx_i] < v_k[idx_k]) {
                idx_i++;
            } else if (v_k[idx_k] < v_i[idx_i]) {
                idx_k++;
            } else {
                u++;
                idx_i++;
                idx_k++;
            }
        }

        int v = 2 * ncol - u;  // count union number

        if(u>0) {
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

// Helper function to project point p onto line segment AB
NumericVector project_point_to_line_segment_cpp(NumericVector p, NumericVector A, NumericVector B) {
  NumericVector AB = B - A;
  double AB_squared = sum(pow(AB, 2.0));
  NumericVector q;

  if (AB_squared == 0.0) {
    q = A;
  } else {
    NumericVector Ap = p - A;
    double t = sum(Ap * AB) / AB_squared;

    if (t < 0.0) {
      q = A;
    } else if (t > 1.0) {
      q = B;
    } else {
      q = A + t * AB;
    }
  }
  return q;
}

// Helper function to project point onto infinite line defined by A and B
NumericVector projPointOnLine_cpp(NumericVector p, NumericVector A, NumericVector B) {
  NumericVector ap = p - A;
  NumericVector ab = B - A;
  double dot_ap_ab = sum(ap * ab);
  double dot_ab_ab = sum(ab * ab);
  if (dot_ab_ab == 0.0) {
    return A;
  }
  return A + (dot_ap_ab / dot_ab_ab) * ab;
}

// [[Rcpp::export]]
List project_point_to_graph(NumericMatrix X,
                            NumericMatrix NodeCoords,
                            List AdjList,
                            IntegerVector ClosestVertex,
                            LogicalVector TipLeaves,
                            bool OrthoProjTip) {

  int N = X.ncol(); // Cells
  int D = X.nrow(); // Dimensions
  NumericMatrix P(D, N);
  IntegerMatrix NearestEdges(N, 2);

  // ClosestVertex is 1-based (from R)
  // AdjList contains 1-based indices (from R)

  for (int i = 0; i < N; ++i) {
    int node_idx_1based = ClosestVertex[i];
    int node_idx = node_idx_1based - 1; // 0-based for C++ indexing into vectors/matrices

    // Check bounds
    if (node_idx < 0 || node_idx >= NodeCoords.ncol()) {
        continue; // Should not happen
    }

    IntegerVector neighbors_1based = AdjList[node_idx];

    double min_dist = R_PosInf;
    NumericVector best_proj(D);
    int best_neighbor = -1;
    bool found_valid = false;

    NumericVector p = X( _, i);
    NumericVector node_coord = NodeCoords( _, node_idx);

    for (int k = 0; k < neighbors_1based.size(); ++k) {
        int neighbor_idx_1based = neighbors_1based[k];
        int neighbor_idx = neighbor_idx_1based - 1;

        if (neighbor_idx == node_idx) continue;
        if (neighbor_idx < 0 || neighbor_idx >= NodeCoords.ncol()) continue;

        NumericVector neighbor_coord = NodeCoords( _, neighbor_idx);
        NumericVector tmp_proj;

        bool is_tip = TipLeaves[node_idx];

        if (is_tip) {
             if (OrthoProjTip) {
                 tmp_proj = projPointOnLine_cpp(p, node_coord, neighbor_coord);
             } else {
                 tmp_proj = project_point_to_line_segment_cpp(p, node_coord, neighbor_coord);
             }
        } else {
             tmp_proj = project_point_to_line_segment_cpp(p, node_coord, neighbor_coord);
        }

        bool valid = true;
        for (int d = 0; d < D; ++d) {
          if (!std::isfinite(tmp_proj[d])) {
            valid = false;
            break;
          }
        }
        if (!valid) {
          tmp_proj = neighbor_coord;
        }

        double dist = sqrt(sum(pow(p - tmp_proj, 2.0)));

        if (dist < min_dist) {
            min_dist = dist;
            best_proj = tmp_proj;
            best_neighbor = neighbor_idx_1based;
            found_valid = true;
        }
    }

    if (found_valid) {
        P( _, i) = best_proj;
        NearestEdges(i, 0) = node_idx_1based;
        NearestEdges(i, 1) = best_neighbor;
    } else {
        // Fallback if no neighbors (isolated node?) - just project to the node itself
        P( _, i) = node_coord;
        NearestEdges(i, 0) = node_idx_1based;
        NearestEdges(i, 1) = node_idx_1based; // Self-loop as fallback
    }
  }

  return List::create(Named("P") = P, Named("nearest_edges") = NearestEdges);
}

// [[Rcpp::export]]
NumericMatrix calc_specificity_cpp(NumericMatrix agg_expr_matrix) {
  int genes = agg_expr_matrix.nrow();
  int groups = agg_expr_matrix.ncol();
  NumericMatrix specificity_mat(genes, groups);

  for (int i = 0; i < genes; ++i) {
    // 1. Match makeprobsvec: p <- p/sum(p); p[is.na(p)] <- 0
    std::vector<double> p(groups);
    std::vector<double> p_raw(groups);
    double sum_raw = 0.0;
    bool has_na = false;
    for (int j = 0; j < groups; ++j) {
      double val = agg_expr_matrix(i, j);
      p_raw[j] = val;
      if (R_IsNA(val) || R_IsNaN(val)) {
        has_na = true;
      } else if (!has_na) {
        sum_raw += val;
      }
    }
    if (has_na) {
      sum_raw = NA_REAL;
    }
    for (int j = 0; j < groups; ++j) {
      double phat = p_raw[j] / sum_raw;
      if (R_IsNA(phat) || R_IsNaN(phat)) {
        phat = 0.0;
      }
      p[j] = phat;
    }

    auto shannon_entropy = [&](const std::vector<double>& v) -> double {
      double min_v = R_PosInf;
      double sum_v = 0.0;
      for (double x : v) {
        if (x < min_v) min_v = x;
        sum_v += x;
      }
      if (min_v < 0.0 || sum_v <= 0.0) {
        return R_PosInf;
      }
      double ent = 0.0;
      for (double x : v) {
        if (x > 0.0) {
          double pn = x / sum_v;
          ent -= pn * std::log2(pn);
        }
      }
      return ent;
    };

    // Check if we can use the optimized path (valid probability vector).
    bool all_finite = true;
    bool any_negative = false;
    double sum_p = 0.0;
    for (double x : p) {
      if (!std::isfinite(x)) all_finite = false;
      if (x < 0.0) any_negative = true;
      sum_p += x;
    }

    if (all_finite && !any_negative && sum_p > 0.0) {
      // Normalize to reduce numerical drift; this matches shannon.entropy's renorm.
      for (int j = 0; j < groups; ++j) {
        p[j] /= sum_p;
      }

      // 2. Calculate Entropy of p
      double H_p = 0.0;
      for (int j = 0; j < groups; ++j) {
        if (p[j] > 0.0) {
          H_p -= p[j] * std::log2(p[j]);
        }
      }

      // 3. Calculate 1 - JS(p, q_j) for each j using optimized formula
      double sum_all_term_k = 0.5 * (-H_p - 1.0);

      for (int j = 0; j < groups; ++j) {
        double term_j_orig = 0.0;
        if (p[j] > 0.0) {
          term_j_orig = 0.5 * (p[j] * std::log2(p[j]) - p[j]);
        }

        double sum_k_neq_j = sum_all_term_k - term_j_orig;

        double m_j = (p[j] + 1.0) / 2.0;
        double term_m_j = 0.0;
        if (m_j > 0.0) {
          term_m_j = m_j * std::log2(m_j);
        }

        double H_m = - (sum_k_neq_j + term_m_j);

        double inner = H_m - 0.5 * H_p; // H(q_j) is 0
        if (inner < 0.0) inner = 0.0;

        double js = std::sqrt(inner);
        specificity_mat(i, j) = 1.0 - js;
      }
    } else {
      // Fallback to exact JS calculation to match R edge cases.
      double H_p = shannon_entropy(p);
      std::vector<double> m(groups);
      for (int j = 0; j < groups; ++j) {
        for (int k = 0; k < groups; ++k) {
          m[k] = (p[k] + (k == j ? 1.0 : 0.0)) / 2.0;
        }
        double H_m = shannon_entropy(m);
        double jsdiv = H_m - 0.5 * H_p; // H(q_j) = 0
        if (std::isinf(jsdiv)) jsdiv = 1.0;
        if (jsdiv < 0.0) jsdiv = 0.0;
        double js = std::sqrt(jsdiv);
        specificity_mat(i, j) = 1.0 - js;
      }
    }
  }
  return specificity_mat;
}
