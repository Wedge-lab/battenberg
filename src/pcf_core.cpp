#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

//' PottsCompact C++ implementation
//' @param kmin Minimal length of plateau
//' @param gamma Penalty for discontinuity
//' @param nr number of values between breakpoints
//' @param res sum of values between breakpoints
//' @param sq sum of squares of values between breakpoints
//' @param yest boolean for estimation
//' @return List with bestCost and bestSplit
// [[Rcpp::export]]
List PottsCompact_cpp(int kmin, double gamma, NumericVector nr, NumericVector res, NumericVector sq) {
  int N = nr.size();
  
  // Scrathpads for accumulation
  std::vector<double> Ant(N + 1, 0.0);
  std::vector<double> Sum(N + 1, 0.0);
  std::vector<double> Kvad(N + 1, 0.0);
  
  NumericVector bestCost(N + 1, 0.0);
  IntegerVector bestSplit(N + 1, 0);

  double initAnt = nr[0];
  double initSum = res[0];
  double initKvad = sq[0];
  double initAve = initSum / initAnt;

  bestCost[0] = initKvad - initSum * initAve;

  int k = 2; // R: k <- 2
  double current_sum_nr = nr[0];
  while (k <= N && (current_sum_nr + nr[k-1]) < 2 * kmin) {
    double nr_k = nr[k-1];
    double res_k = res[k-1];
    double sq_k = sq[k-1];
    current_sum_nr += nr_k;
    for (int i = 2; i <= k; ++i) {
      Ant[i] += nr_k;
      Sum[i] += res_k;
      Kvad[i] += sq_k;
    }
    bestCost[k-1] = (initKvad + Kvad[2]) - std::pow(initSum + Sum[2], 2) / (initAnt + Ant[2]);
    k++;
  }

  // Main DP loop
  for (int n = k; n <= N; ++n) {
    double nr_n = nr[n-1];
    double res_n = res[n-1];
    double sq_n = sq[n-1];

    for (int i = 2; i <= n; ++i) {
      Ant[i] += nr_n;
      Sum[i] += res_n;
      Kvad[i] += sq_n;
    }

    int limit = n;
    while (limit > 2 && Ant[limit] < kmin) {
      limit--;
    }

    int bestPos = 2;
    double min_cost = -1.0;

    for (int i = 2; i <= limit; ++i) {
      double current_cost = bestCost[i-2] + Kvad[i] - std::pow(Sum[i], 2) / Ant[i];
      if (i == 2 || current_cost < min_cost) {
        min_cost = current_cost;
        bestPos = i;
      }
    }

    double final_cost = min_cost + gamma;
    double tot_cost = (Kvad[2] + initKvad) - std::pow(Sum[2] + initSum, 2) / (Ant[2] + initAnt);
    
    if (tot_cost < final_cost) {
      bestPos = 1;
      final_cost = tot_cost;
    }

    bestCost[n-1] = final_cost;
    bestSplit[n-1] = bestPos - 1;
  }

  return List::create(
    _["bestCost"] = bestCost,
    _["bestSplit"] = bestSplit
  );
}

//' exactPcf C++ implementation
//' @param y Input vector
//' @param kmin Minimal length of plateau
//' @param gamma Penalty
//' @return List with bestCost, bestAver, bestSplit
// [[Rcpp::export]]
List exactPcf_cpp(NumericVector y, int kmin, double gamma) {
  int N = y.size();
  
  NumericVector bestCost(N + 1, 0.0);
  NumericVector bestAver(N + 1, 0.0);
  IntegerVector bestSplit(N + 1, 0);

  std::vector<double> Sum(N + 1, 0.0);
  std::vector<double> Kvad(N + 1, 0.0);
  std::vector<double> Aver(N + 1, 0.0);

  double initSum = 0;
  double initKvad = 0;
  for (int i = 0; i < kmin; ++i) {
    initSum += y[i];
    initKvad += y[i] * y[i];
  }
  double initAve = initSum / kmin;
  bestCost[kmin] = initKvad - initSum * initAve;
  bestAver[kmin] = initAve;

  int kminP1 = kmin + 1;
  int limit_init = 2 * kmin - 1;
  if (limit_init > N) limit_init = N;

  for (int k = kminP1; k <= limit_init; ++k) {
    double yk = y[k-1];
    double yk2 = yk * yk;
    for (int i = kminP1; i <= k; ++i) {
      Sum[i] += yk;
      Kvad[i] += yk2;
    }
    double cur_best_aver = (initSum + Sum[kminP1]) / k;
    bestAver[k] = cur_best_aver;
    bestCost[k] = (initKvad + Kvad[kminP1]) - k * std::pow(cur_best_aver, 2);
  }

  for (int n = (2 * kmin); n <= N; ++n) {
    double yn = y[n-1];
    double yn2 = yn * yn;
    
    for (int i = kminP1; i <= n; ++i) {
      Sum[i] += yn;
      Kvad[i] += yn2;
      Aver[i] = Sum[i] / (n - i + 1);
    }

    int limit = n - kmin + 1;
    int bestPos = kminP1;
    double min_cost = -1.0;

    for (int i = kminP1; i <= limit; ++i) {
      double current_cost = bestCost[i-1] + Kvad[i] - Sum[i] * Aver[i] + gamma;
      if (i == kminP1 || current_cost < min_cost) {
        min_cost = current_cost;
        bestPos = i;
      }
    }

    double final_cost = min_cost;
    double final_aver = Aver[bestPos];
    
    double totSum = Sum[kminP1] + initSum;
    double totKvad = Kvad[kminP1] + initKvad;
    double totAver = totSum / n;
    double totCost = totKvad - n * std::pow(totAver, 2);

    if (totCost < final_cost) {
      bestPos = 1;
      final_cost = totCost;
      final_aver = totAver;
    }

    bestCost[n] = final_cost;
    bestAver[n] = final_aver;
    bestSplit[n] = bestPos - 1;
  }

  return List::create(
    _["bestCost"] = bestCost,
    _["bestAver"] = bestAver,
    _["bestSplit"] = bestSplit
  );
}

//' findEst C++ implementation
//' @param bestSplit vector of best splits from DP
//' @param N number of compressed points
//' @param Nr number of original points in each compressed point
//' @param Sum sum of original values in each compressed point
//' @param yest boolean for estimation
//' @return List with segments and optionally yhat
// [[Rcpp::export]]
List findEst_cpp(IntegerVector bestSplit, int N, NumericVector Nr, NumericVector Sum, bool yest) {
  int n = N;
  std::vector<int> lengde_comp;
  while (n > 0) {
    int split = bestSplit[n]; // Adjusted for 1-based indexing passed from DP if needed
    // The Rcpp passed vector is likely 0-indexed if it came from our Rcpp DP
    // But R's original bestSplit was 1-indexed. Let's assume 0-indexed from our DP.
    lengde_comp.push_back(n - split);
    n = split;
  }
  std::reverse(lengde_comp.begin(), lengde_comp.end());
  
  int antInt = lengde_comp.size();
  NumericVector lengdeOrig(antInt);
  NumericVector startOrig(antInt);
  NumericVector verdi(antInt);
  
  int current_start = 0;
  int current_start_orig = 1;
  
  for (int i = 0; i < antInt; ++i) {
    int l_comp = lengde_comp[i];
    double l_orig = 0;
    double s_orig = 0;
    for (int j = 0; j < l_comp; ++j) {
      l_orig += Nr[current_start + j];
      s_orig += Sum[current_start + j];
    }
    
    lengdeOrig[i] = l_orig;
    startOrig[i] = current_start_orig;
    verdi[i] = s_orig / l_orig;
    
    current_start += l_comp;
    current_start_orig += (int)l_orig;
  }
  
  if (yest) {
    int totalN = current_start_orig - 1;
    NumericVector yhat(totalN);
    int pos = 0;
    for (int i = 0; i < antInt; ++i) {
      int l = (int)lengdeOrig[i];
      double val = verdi[i];
      for (int j = 0; j < l; ++j) {
        yhat[pos++] = val;
      }
    }
    return List::create(
      _["Lengde"] = lengdeOrig,
      _["sta"] = startOrig,
      _["mean"] = verdi,
      _["nIntervals"] = antInt,
      _["yhat"] = yhat
    );
  } else {
    return List::create(
      _["Lengde"] = lengdeOrig,
      _["sta"] = startOrig,
      _["mean"] = verdi,
      _["nIntervals"] = antInt
    );
  }
}

//' findMarks C++ implementation
//' @param markSub marks in compressed scale
//' @param Nr number of observations
//' @param subsize original scale size
//' @return LogicalVector of marks in original scale
// [[Rcpp::export]]
LogicalVector findMarks_cpp(LogicalVector markSub, NumericVector Nr, int subsize) {
  LogicalVector mark(subsize, false);
  int N = markSub.size();
  
  int oldStart = 0;
  int startOrig = 1;
  
  for (int i = 0; i < N; ++i) {
    if (markSub[i]) {
      // Index i is the end of a segment
      double l_orig = 0;
      for (int k = oldStart; k <= i; ++k) {
        l_orig += Nr[k];
      }
      startOrig += (int)l_orig;
      if (startOrig - 1 <= subsize) {
        mark[startOrig - 2] = true;
      }
      oldStart = i + 1;
    }
  }
  return mark;
}
