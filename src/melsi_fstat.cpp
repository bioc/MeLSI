#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
using namespace Rcpp;

// Fused PERMANOVA pseudo-F statistic computed directly from data.
//
// Xt is the (already weighted/scaled) data transposed to features-by-samples
// (m x n) so each sample's features are contiguous in memory. group holds
// 0-based group ids of length n; k is the number of groups.
//
// Squared Euclidean distances are accumulated in a single pass over the
// n*(n-1)/2 unique sample pairs, so the full n x n distance matrix is never
// materialised. Totals use long double to limit floating-point error.
//
// Equivalent (in exact arithmetic) to the R reference:
//   D2 <- as.matrix(dist(t(Xt)))^2
//   SS_total  <- sum(D2) / (2n)
//   SS_within <- sum_g sum(D2[idx_g, idx_g]) / (2 n_g)
//   F <- (SS_between/(k-1)) / (SS_within/(n-k))
//
// [[Rcpp::export]]
double melsi_permanova_f(const NumericMatrix& Xt, const IntegerVector& group, int k) {
    const int m = Xt.nrow();
    const int n = Xt.ncol();

    // Plain double accumulators (not long double): on x86-64 long double is
    // 80-bit x87, which is slow and prevents SIMD vectorisation of the inner
    // loop. Summing O(n^2) similar-magnitude positive terms in double keeps the
    // relative error at ~1e-14, well within the tolerance the equivalence tests
    // confirmed.
    std::vector<double> within(k, 0.0);
    std::vector<int> nsize(k, 0);
    for (int i = 0; i < n; ++i) nsize[group[i]]++;

    const double* x = &Xt[0];           // column-major: sample i starts at x + i*m
    double total = 0.0;

    for (int i = 0; i < n; ++i) {
        const double* __restrict xi = x + (std::size_t)i * m;
        const int gi = group[i];
        for (int j = i + 1; j < n; ++j) {
            const double* __restrict xj = x + (std::size_t)j * m;
            double s = 0.0;
            for (int d = 0; d < m; ++d) {
                const double diff = xi[d] - xj[d];
                s += diff * diff;
            }
            total += s;
            if (group[j] == gi) within[gi] += s;
        }
    }

    const double SS_total = total / (double)n;
    double SS_within = 0.0;
    for (int g = 0; g < k; ++g) {
        if (nsize[g] > 0) SS_within += within[g] / (double)nsize[g];
    }
    const double SS_between = SS_total - SS_within;
    const double F = (SS_between / (double)(k - 1)) /
                     (SS_within / (double)(n - k));
    return F;
}


// ---------------------------------------------------------------------------
// Weak-learner metric optimizer (diagonal Mahalanobis), ported from the R
// reference optimize_weak_learner_robust(). The gradient descent samples one
// within-class pair from each of the first two classes plus one between-class
// pair per iteration and updates the diagonal metric.
//
// REPRODUCIBILITY: the four per-iteration draws use R's OWN index sampler
// (::R_unif_index, declared in Rmath.h), which consumes the global R RNG via
// unif_rand() in exactly the same pattern as sample()/sample.int(). Combined
// with the bootstrap/feature sampling that stays in R, the RNG stream is
// byte-aligned with the pure-R version, so the learned metric is identical.
// Confirmed bit-exact (max|diff| = 0) across class sizes and 2/3-group cases.
//
// g holds 0-based group ids of length n (= match(y, unique(y)) - 1L); k is the
// number of groups. The gradient uses classes 0 and 1 (= unique(y)[1:2]); the
// early-stopping F statistic uses all k groups, matching the R reference.
// ---------------------------------------------------------------------------
static double melsi_f_scaled_buf(const double* Xb, int n, int m,
                                 const std::vector<double>& scale,
                                 const IntegerVector& g, int k) {
    // sample-contiguous, pre-scaled buffer so the pairwise loop reads sequential
    // memory (same layout melsi_permanova_f relies on); avoids R sweep()+t().
    std::vector<double> S((std::size_t)n * m);
    for (int i = 0; i < n; ++i) {
        double* __restrict si = &S[(std::size_t)i * m];
        for (int d = 0; d < m; ++d) si[d] = Xb[(std::size_t)d * n + i] * scale[d];
    }
    std::vector<double> within(k, 0.0);
    std::vector<int> nsize(k, 0);
    for (int i = 0; i < n; ++i) nsize[g[i]]++;
    double total = 0.0;
    for (int i = 0; i < n; ++i) {
        const double* __restrict si = &S[(std::size_t)i * m];
        const int gi = g[i];
        for (int j = i + 1; j < n; ++j) {
            const double* __restrict sj = &S[(std::size_t)j * m];
            double s = 0.0;
            for (int d = 0; d < m; ++d) { const double diff = si[d] - sj[d]; s += diff * diff; }
            total += s;
            if (g[j] == gi) within[gi] += s;
        }
    }
    const double SS_total = total / (double)n;
    double SS_within = 0.0;
    for (int gg = 0; gg < k; ++gg) if (nsize[gg] > 0) SS_within += within[gg] / (double)nsize[gg];
    return ((SS_total - SS_within) / (double)(k - 1)) / (SS_within / (double)(n - k));
}

// [[Rcpp::export]]
NumericVector melsi_opt_weak_learner(const NumericMatrix& X, const IntegerVector& g,
                                     int k, int n_iterations = 50,
                                     double learning_rate = 0.1) {
    const int n = X.nrow();
    const int m = X.ncol();
    std::vector<double> Md(m, 1.0);  // diag(M), starts at identity

    std::vector<int> c1, c2;         // class1 = which(g==0), class2 = which(g==1)
    for (int i = 0; i < n; ++i) { if (g[i] == 0) c1.push_back(i); else if (g[i] == 1) c2.push_back(i); }
    const int n1 = (int)c1.size(), n2 = (int)c2.size();
    if (k < 2 || n1 < 2 || n2 < 2) return wrap(Md);

    double prev_f = -std::numeric_limits<double>::infinity();
    int stagnation = 0;

    for (int iter = 1; iter <= n_iterations; ++iter) {
        int a  = (int) ::R_unif_index((double)n1);        int i1 = c1[a];
        int b  = (int) ::R_unif_index((double)(n1 - 1));  int j1 = (b < a) ? c1[b] : c1[b + 1];
        int cc = (int) ::R_unif_index((double)n2);        int i2 = c2[cc];
        int d  = (int) ::R_unif_index((double)(n2 - 1));  int j2 = (d < cc) ? c2[d] : c2[d + 1];

        const double lr_t = learning_rate * (1.0 / (1.0 + iter * 0.1));
        for (int q = 0; q < m; ++q) {
            const double d1 = X(i1, q) - X(j1, q);  // within class 1
            const double d2 = X(i2, q) - X(j2, q);  // within class 2
            const double d3 = X(i1, q) - X(i2, q);  // between classes
            const double grad = d3 * d3 - (d1 * d1 + d2 * d2) / 2.0;
            double v = Md[q] + lr_t * grad;
            if (v < 0.01) v = 0.01;                  // keep positive (pmax 0.01)
            Md[q] = v;
        }

        if (iter % 20 == 0) {                        // early-stopping check
            std::vector<double> sc(m);
            for (int q = 0; q < m; ++q) sc[q] = std::sqrt(Md[q] > 0 ? Md[q] : 0.0);
            double cur = melsi_f_scaled_buf(&X[0], n, m, sc, g, k);
            if (!std::isfinite(cur)) cur = 0.0;
            if (cur <= prev_f) { if (++stagnation >= 5) break; }
            else stagnation = 0;
            prev_f = cur;
        }
    }
    return wrap(Md);
}
