//=========================================================================================
#include <Rcpp.h>
using namespace Rcpp;
//=========================================================================================
#define tiny 1e-30
#define epsilon 1e-6
// functions used in smpqrs
NumericMatrix bsplines(NumericVector x, NumericVector t, int degree, int derivative);
NumericVector bspline(NumericVector x, NumericVector t, int degree, int j, int derivative);
double bsp(double x, NumericVector t, int degree, int j);
double dbsp(double x, NumericVector t, int degree, int j, int derivative);
NumericMatrix jump_bsplines(NumericVector t, int degree);
NumericVector dim2knots(NumericVector predictor, int dimension, int degree);
NumericVector knots2t(NumericVector knots, int degree);
NumericVector lambdas_all_alpha(int number_lambdas_alpha,
                                double lambda_alpha_max, double epsilon_lambda);
NumericVector lambdas_all_beta(int number_lambdas,
                               double lambda_beta_max, double epsilon_lambda);
// array operations
ListOf<NumericVector> NumericVectors(int size);
ListOf<IntegerVector> IntegerVectors(int size);
IntegerVector support_of_vector(NumericVector v);
IntegerVector bubble_order(NumericVector vec);
NumericVector subsetNumVec(NumericVector x, IntegerVector index);
IntegerVector subsetIntVec(IntegerVector x, IntegerVector index);
NumericVector knot_candidate(NumericMatrix predictors, int knot_size);
double norm_2(NumericVector u);
double norm_square(NumericVector u);
double zlambda(NumericVector a, double b, double c, NumericVector d, double lambda);
NumericVector delta(NumericVector d);
double q_lambda_point(double z, NumericVector a, double b, double c, NumericVector d, double lambda);
NumericVector alpha2xi(NumericVector alpha);
NumericVector xi2alpha(NumericVector xi);
NumericVector deriv_xi_j(NumericVector xi, NumericVector alpha, int j);
double atan4(double y, double x);
//=========================================================================================
// Main function
// Single Index Model using B-splines and Polar Coordinate
// [[Rcpp::export]]
List BPSI(NumericVector responses,
          NumericMatrix predictors,
          NumericVector initial_alpha,
          int           degree,
          int           number_interior_knots,
          int           number_lambdas_alpha,
          int           number_lambdas_beta,
          double        lambda_alpha_max = 100,
          double        lambda_beta_max = 10,
          double        epsilon_lambda = 1e-6,
          int           maxiter = 200,
          double        epsilon_iterations = 1e-4)
{
   Rcout << "=============================================================\n";
   Rcout << "Single Index Model using B-splines and Spherical Coordinates\n";
   Rcout << "Version 1.0 by SDMLAB (January 17, 2018)\n";
   Rcout << "Department of Statistics, Korea University, Korea\n";
   Rcout << "=============================================================\n";
   List results(number_lambdas_alpha + 1);
   List results_beta(number_lambdas_beta);
   // setup
   int sample_size = responses.size();
   int order = degree + 1;
   // initial knots
   NumericVector initial_knots = knot_candidate(predictors, number_interior_knots);
   NumericVector initial_t = knots2t(initial_knots, degree);
   int initial_dimension = initial_t.size() - degree - 1;
   int dimension = initial_dimension;
   int number_penalty = dimension - order;
   int pm1 = initial_alpha.size();
   int p = pm1 + 1;
   int active_pm1 = initial_alpha.size();
   int active_p = active_pm1 + 1;
   int j, k, m;
   int iter, iter_alpha, iter_beta;
   IntegerVector xi_index(p);
   //IntegerVector active_xi_index(p);
   // initial coefficients
   NumericVector initial_xi = alpha2xi(initial_alpha);
   NumericVector xi = clone(initial_xi);
   NumericVector active_xi = clone(initial_xi);
   NumericVector active_alpha = clone(initial_alpha);
   // initial dot_products
   NumericVector initial_dot_products(sample_size);
   for(m = 0; m < p; m++)
      initial_dot_products += initial_xi[m] * predictors(_, m);
   NumericMatrix initial_basis = bsplines(initial_dot_products, initial_t, degree, 0);
   NumericMatrix initial_jump = jump_bsplines(initial_t, degree);
   NumericVector residuals = clone(responses);
   NumericVector fitted_value(sample_size);
   NumericVector bj(sample_size);
   NumericVector am(sample_size);
   NumericVector partial_residuals(sample_size);
   NumericVector delta(sample_size);
   double b, c, R, R_beta, lambda_alpha, lambda_beta, solution_alpha, alpha_star_points, alpha_star;
   double store_R = R_PosInf;
   double store_R_beta = R_PosInf;
   double store_store_R_beta = R_PosInf;
   double alpha_star_candidate_points;
   // penalty alpha
   NumericVector candidate_alpha(1);
   NumericVector weight_alpha(1);
   weight_alpha[0] = 1;
   NumericVector alpha_penalty_constants(pm1);
   // all lambdas
   NumericVector lambdas_alpha = lambdas_all_alpha(number_lambdas_alpha,
                                                   lambda_alpha_max, epsilon_lambda);
   NumericVector lambdas_beta = lambdas_all_beta(number_lambdas_beta,
                                                 lambda_beta_max, epsilon_lambda);
   // aic, bic
   NumericMatrix bic_matrix(number_lambdas_alpha, number_lambdas_beta);
   NumericMatrix aic_matrix(number_lambdas_alpha, number_lambdas_beta);
   // Module FIT
   for (int lambda_index_alpha = 0; lambda_index_alpha < number_lambdas_alpha; lambda_index_alpha++)
   {
      lambda_alpha = lambdas_alpha[lambda_index_alpha];
      // initialized : before pruning !!!
      NumericVector knots = clone(initial_knots);
      NumericVector t = clone(initial_t);
      dimension = initial_dimension;
      number_penalty = dimension - order;
      // initial coefficients
      NumericVector alpha = clone(initial_alpha);
      NumericVector beta(dimension);
      NumericVector xi = clone(initial_xi);
      NumericVector active_xi = clone(initial_xi);
      p = xi.size();
      pm1 = p - 1;
      active_p = xi.size();
      active_pm1 = active_p - 1;
      IntegerVector xi_index = seq(0, p - 1);
      IntegerVector active_xi_index = seq(0, p - 1);
      // initial dot_products
      NumericVector dot_products = clone(initial_dot_products);
      NumericMatrix basis = clone(initial_basis);
      NumericMatrix jump = clone(initial_jump);
      residuals = clone(responses);
      // fitted_value = clone(responses);
      R = 0.0;
      store_R = R_PosInf;
      store_R_beta = R_PosInf;
      store_store_R_beta = R_PosInf;
      // lambda_alpha
      for (int lambda_index_beta = 0; lambda_index_beta < number_lambdas_beta; lambda_index_beta++)
      {
         Rcout << "\r" << "lambda index (alpha, beta) : (" <<
            lambda_index_alpha << "," << lambda_index_beta << ") / (" <<
               number_lambdas_alpha << "," << number_lambdas_beta << ")";
         lambda_beta = lambdas_beta[lambda_index_beta];
         // iteration 0. until converge
         for (iter = 0; iter < 5; iter++)
         {
            // iteration 1. until beta converge
            for (iter_beta = 0; iter_beta < maxiter; iter_beta++)
            {
               // Module UPDATE
               // Module UPDATE beta
               for (j = 0; j < dimension; j++)
               {
                  //a = bspline(dot_products, t, degree, j, 0);
                  bj = basis(_, j);
                  partial_residuals = residuals + beta[j] * bj;
                  b = sum(bj * bj);
                  c = sum(bj * partial_residuals) / (b + tiny);
                  if (dimension == order)
                     beta[j] = c;
                  else
                  {
                     // compute a, d vector for penalty term
                     // d_1 |beta_j - a_1| + ... + d_M |beta_j - a_M|
                     NumericVector rowjump_j = jump(j, _);
                     NumericVector d = rowjump_j[abs(rowjump_j) > epsilon];
                     int nonzero_size = d.size();
                     if (nonzero_size == 0)
                        beta[j] = c;
                     else
                     {
                        // compute a
                        NumericVector a(nonzero_size);
                        IntegerVector nonzero_index = seq(0, rowjump_j.size() - 1);
                        nonzero_index = nonzero_index[abs(rowjump_j) > epsilon];
                        for (k = 0; k < nonzero_size; k++)
                           a[k] = beta[j] - sum(jump(_, nonzero_index[k]) * beta) / d[k];
                        // abs d after compute a
                        d = abs(d);
                        // update beta_j
                        beta[j] = zlambda(a, b, c, d, lambda_beta);
                     }
                  }
                  // update residuals
                  residuals = partial_residuals - beta[j] * bj;
               }
               // Module Prune
               if (number_penalty > 0)
               {
                  NumericVector penalty(number_penalty);
                  for (k = 0; k < number_penalty; k++)
                     penalty[k] = sum(jump(_, k) * beta);
                  IntegerVector penalty_check(number_penalty);
                  penalty_check[abs(penalty) < epsilon] = 1;
                  if (sum(penalty_check) > 0)
                  {
                     // re-compute
                     IntegerVector prune_coef(dimension);
                     IntegerVector prune_knots(number_penalty + 2);
                     IntegerVector prune_index = seq(0, number_penalty - 1);
                     prune_index = prune_index[abs(penalty) < epsilon];
                     prune_coef[prune_index] = 1;
                     prune_knots[prune_index + 1] = 1;
                     beta = beta[prune_coef == 0];
                     knots = knots[prune_knots == 0];
                     dimension = beta.size();
                     // fitted_values
                     number_penalty = dimension - order;
                     t = knots2t(knots, degree);
                     basis = bsplines(dot_products, t, degree, 0);
                     if (number_penalty > 0)
                        jump = jump_bsplines(t, degree);
                  }
               }
               // UPDATE fitted_value
               fitted_value = rep(0, sample_size);
               for (j = 0; j < dimension; j ++)
                  fitted_value += beta[j] * basis(_, j);
               residuals = responses - fitted_value;
               // Module CHECK_CONVERGENCE
               R_beta = 0.5 * sum(pow(residuals, 2.0));
               for (k = 0; k < number_penalty; k++)
                  R_beta += lambda_beta * std::abs(sum(jump(_, k) * beta));
               if (std::abs(R_beta - store_R_beta) < epsilon_iterations)
                  break;
               // save R
               store_R_beta = R_beta;
            }
            if (std::abs(R_beta - store_store_R_beta) < epsilon)
               break;
            // save R
            store_store_R_beta = R_beta;
            // iteration 2. until alpha converge
            for (iter_alpha = 0; iter_alpha < 1; iter_alpha++)
            {
               if (active_p == 1)
                  break;
               NumericMatrix basis_p = bsplines(dot_products, t, degree, 1);
               // Module UPDATE alpha
               // alpha_1 ~ alpha_pm1 - 1
               for (m = 0; m < active_pm1 - 1; m++)
               {
                  am = rep(0, sample_size);
                  delta = rep(0, sample_size);
                  NumericVector deriv_xi = deriv_xi_j(active_xi, alpha[seq(0, active_pm1 - 1)], m);
                  for (k = 0; k < active_p; k++)
                  {
                     j = active_xi_index[k];
                     delta += predictors(_, j) * deriv_xi[k];
                  }
                  for (j = 0; j < dimension; j++)
                     am += beta[j] * basis_p(_, j) * delta;
                  partial_residuals = residuals + alpha[m] * am;
                  b = sum(am * am);
                  c = sum(am * partial_residuals) / (b + tiny);
                  alpha[m] = c;
                  // penalty alpha : localized penalty
                  alpha_star_points = R_PosInf;
                  for (int i = 0; i < 3; i++)
                  {
                     candidate_alpha[0] = PI * i / 2;
                     solution_alpha = zlambda(candidate_alpha, b, c, weight_alpha, lambda_alpha);
                     alpha_star_candidate_points = q_lambda_point(solution_alpha, candidate_alpha, b, c, weight_alpha, lambda_alpha);
                     if (alpha_star_candidate_points < alpha_star_points)
                     {
                        alpha_star_points = alpha_star_candidate_points;
                        alpha_star = solution_alpha;
                        alpha_penalty_constants[m] = candidate_alpha[0];
                     }
                  }
                  alpha[m] = alpha_star;
                  // UPDATE residuals
                  residuals = partial_residuals - alpha[m] * am;
               }
               // alpha_pm1
               am.fill(0.0);
               NumericVector deriv_xi = deriv_xi_j(active_xi, alpha[seq(0, active_pm1 - 1)], active_pm1 - 1);
               delta.fill(0.0);
               for (k = 0; k < active_p; k++)
               {
                  j = active_xi_index[k];
                  delta += predictors(_, j) * deriv_xi[k];
               }
               for (j = 0; j < dimension; j++)
                  am += beta[j] * basis_p(_, j) * delta;
               partial_residuals = residuals + alpha[active_pm1 - 1] * am;
               b = sum(am * am);
               c = sum(am * partial_residuals) / (b + tiny);
               alpha[active_pm1 - 1] = alpha_star;
               // penalty alpha : localized penalty
               alpha_star_points = R_PosInf;
               for (int i = 0; i < 5; i++)
               {
                  candidate_alpha[0] = PI * i / 2;
                  solution_alpha = zlambda(candidate_alpha, b, c, weight_alpha, lambda_alpha);
                  alpha_star_candidate_points = q_lambda_point(solution_alpha, candidate_alpha, b, c, weight_alpha, lambda_alpha);
                  if (alpha_star_candidate_points < alpha_star_points)
                  {
                     alpha_star_points = alpha_star_candidate_points;
                     alpha_star = solution_alpha;
                     alpha_penalty_constants[active_pm1 - 1] = candidate_alpha[0];
                  }
               }
               alpha[active_pm1 - 1] = alpha_star;
               // UPDATE xi
               xi.fill(0.0);
               xi[active_xi_index] = alpha2xi(alpha[seq(0, active_pm1 - 1)]);
               // Module Prune_alpha
               xi_index = seq(0, p - 1);
               active_xi_index = xi_index[abs(xi) > tiny];
               active_xi = xi[active_xi_index];
               active_p = active_xi.size();
               active_pm1 = active_p - 1;
               // UPDATE dot_products
               dot_products.fill(0.0);
               for(m = 0; m < active_p; m++)
               {
                  k = active_xi_index[m];
                  dot_products += active_xi[m] * predictors(_, k);
               }
               // UPDATE fitted_value
               fitted_value = rep(0, sample_size);
               for (j = 0; j < dimension; j ++)
                  fitted_value += beta[j] * bspline(dot_products, t, degree, j, 0);
               residuals = responses - fitted_value;
               if (active_p > 1)
                  alpha[seq(0, active_pm1 - 1)] = xi2alpha(active_xi);
            }
            basis = bsplines(dot_products, t, degree, 0);
         }
         // Rcout << "iter = " << iter << std::endl;
         // Results
         bic_matrix(lambda_index_alpha, lambda_index_beta) = sample_size * log(sum(pow(residuals, 2.0)) / sample_size) + (dimension + active_p) * log(sample_size);
         aic_matrix(lambda_index_alpha, lambda_index_beta) = sample_size * log(sum(pow(residuals, 2.0)) / sample_size) + (dimension + active_p) * 2.0;
         // final alpha ?
         //alpha = xi2alpha(xi);
         results_beta[lambda_index_beta] = List::create(_["beta"] = clone(beta),
                                                        _["xi"] = clone(xi),
                                                        _["dimension"] = dimension,
                                                        _["knots"] = clone(knots),
                                                        _["t"] = clone(t),
                                                        _["dot_products"] = clone(dot_products),
                                                        _["fitted_value"] = clone(fitted_value),
                                                        _["active_xi_dimension"] = active_p,
                                                        _["active_xi_index"] = clone(active_xi_index) + 1);
      }
      results[lambda_index_alpha] = clone(results_beta);
   }
   results[number_lambdas_alpha] = List::create(_["bic_matrix"] = bic_matrix,
                                                _["aic_matrix"] = aic_matrix,
                                                _["lambdas_alpha"] = lambdas_alpha,
                                                _["lambdas_beta"]    = lambdas_beta);
   Rcout << "\n";
   return results;
}
//=========================================================================================
// [[Rcpp::export]]
NumericMatrix bsplines(NumericVector x, NumericVector t,
                       int degree = 1, int derivative = 0)
{
   int x_length = x.size();
   int dimension = t.size() - degree - 1;
   NumericMatrix bm(x_length, dimension);
   for (int j = 0; j < dimension; j++)
      bm(_, j) = bspline(x, t, degree, j, derivative);
   return bm;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector bspline(NumericVector x, NumericVector t,
                      int degree = 1, int j = 0, int derivative = 0)
{
   int x_length = x.size();
   NumericVector b(x_length);
   int i;
   int k = degree + 1; // k = order = degree + 1
   if (derivative == 0)
   {
      for (i = 0; i < x_length; i++)
      {
         b[i] = 0;
         if ((t[j] <= x[i]) && (x[i] < t[j + k]))
            b[i] = bsp(x[i], t, degree, j);
      }
   }
   else
   {
      for (i = 0; i < x_length; i++)
      {
         b[i] = 0;
         if ((t[j] <= x[i]) && (x[i] < t[j + k]))
            b[i] = dbsp(x[i], t, degree, j, derivative);
      }
   }
   return b;
}
//=========================================================================================
// [[Rcpp::export]]
double bsp(double x, NumericVector t, int degree, int j)
{
   if (degree == 0)
   {
      if ((t[j] <= x) && (x < t[j + 1]))
         return 1;
      else
         return 0;
   }
   else
   {
      double a, b, c, d;
      int k = degree + 1; // k = order = degree + 1
      int jd = j + degree, jpk = j + k, jp1 = j + 1;
      c = t[jd] - t[j];
      if (c > 0)
         a = (x - t[j]) / c;
      else
         a = 0;
      d = t[jpk] - t[jp1];
      if (d > 0)
         b = (t[jpk] - x) / (t[jpk] - t[jp1]);
      else
         b = 0;
      return a * bsp(x, t, degree - 1, j) + b * bsp(x, t, degree - 1, jp1);
   }
}
//=========================================================================================
// [[Rcpp::export]]
double dbsp(double x, NumericVector t, int degree, int j, int derivative)
{
   if (derivative == 0)
      return bsp(x, t, degree, j);
   else
   {
      double a, b, c, d;
      c = t[j + degree] - t[j];
      if (c > 0)
         a = degree / c;
      else
         a = 0;
      int k = degree + 1; // k = order = degree + 1
      int jp1 = j + 1;
      d = t[j + k] - t[jp1];
      if (d > 0)
         b = degree / d;
      else
         b = 0;
      return a * dbsp(x, t, degree - 1, j,     derivative - 1) -
         b * dbsp(x, t, degree - 1, j + 1, derivative - 1);
   }
}
//=========================================================================================
// [[Rcpp::export]]
NumericMatrix jump_bsplines(NumericVector t, int degree = 1)
{
   int k = degree + 1; // k = order = degree + 1
   int dimension = t.size() - k;
   NumericMatrix jump(dimension, dimension - k);
   NumericVector x(dimension - k + 1);
   NumericVector derivative(dimension);
   int j, l;
   for (j = 0; j < x.size(); j++)
      x[j] = 0.5 * (t[j + k - 1] + t[j + k]);
   for (j = 0; j < dimension; j++)
   {
      derivative = bspline(x, t, degree, j, degree);
      for (l = 0; l < (dimension - k); l++)
         jump(j, l) = derivative[l + 1] - derivative[l];
   }
   return jump;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector knots2t(NumericVector knots, int degree = 1)
{
   double d = mean(diff(knots));
   double min_knots = min(knots);
   double max_knots = max(knots);
   int number_knots = knots.size();
   NumericVector t(2 * degree + number_knots);
   int j;
   for (j = 0; j < degree; j++)
      t[j] = min_knots - d * (degree - j);
   int k = j;
   for (j = 0; j < number_knots; j++)
   {
      t[k] = knots[j];
      k++;
   }
   for (j = 0; j < degree; j++)
   {
      t[k] = max_knots + d * (j + 1);
      k++;
   }
   return t;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector dim2knots(NumericVector predictor, int dimension, int degree)
{
   int sample_size = predictor.size();
   NumericVector knots(dimension - degree + 1);
   int knot_size = knots.size();
   int j;
   double i;
   for (j = 0; j < (knot_size - 1); j++)
   {
      i = sample_size / (knot_size - 1.0) * j;
      knots[j] = predictor[i];
   }
   knots[knot_size - 1] = max(predictor) + 1e-5;
   return knots;
}
//=========================================================================================
// [[Rcpp::export]]
ListOf<NumericVector> NumericVectors(int size)
{
   List lov(size);
   NumericVector v;
   for (int i = 0; i < size; i++)
      lov[i] = v;
   return lov;
}
//=========================================================================================
// [[Rcpp::export]]
ListOf<IntegerVector> IntegerVectors(int size)
{
   List lov(size);
   NumericVector v;
   for (int i = 0; i < size; i++)
      lov[i] = v;
   return lov;
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector support_of_vector(NumericVector v)
{
   int n = v.size();
   IntegerVector z2n(n);
   for (int i = 0; i < n; i++)
      z2n[i] = i;
   return z2n[abs(v) > 1e-6];
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector bubble_order(NumericVector vec)
{
   double tmp = 0;
   NumericVector clone_vec = clone(vec);
   int n = vec.size();
   IntegerVector outvec = seq(1, n);
   int itmp, no_swaps, passes;
   passes = 0;
   while(true)
   {
      no_swaps = 0;
      for (int i = 0; i < n - 1 - passes; ++i)
      {
         if(clone_vec[i] > clone_vec[i+1])
         {
            no_swaps++;
            tmp = clone_vec[i];
            clone_vec[i] = clone_vec[i + 1];
            clone_vec[i + 1] = tmp;
            itmp = outvec[i];
            outvec[i] = outvec[i+1];
            outvec[i+1] = itmp;
         }
      }
      if (no_swaps == 0)
         break;
      passes++;
   }
   return outvec;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector subsetNumVec(NumericVector x, IntegerVector index)
{
   // Length of the index vector
   int n = index.size();
   // Initialize output vector
   NumericVector out(n);
   // Subtract 1 from index as C++ starts to count at 0
   index = index - 1;
   // Loop through index vector and extract values of x at the given positions
   for (int i = 0; i < n; i++)
      out[i] = x[index[i]];
   // Return output
   return out;
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector subsetIntVec(IntegerVector x, IntegerVector index)
{
   // Length of the index vector
   int n = index.size();
   // Initialize output vector
   IntegerVector out(n);
   // Subtract 1 from index as C++ starts to count at 0
   index = index - 1;
   // Loop through index vector and extract values of x at the given positions
   for (int i = 0; i < n; i++)
      out[i] = x[index[i]];
   // Return output
   return out;
}
//=========================================================================================
// knot_candidate
//[[Rcpp::export]]
NumericVector knot_candidate(NumericMatrix predictors, int knot_size)
{
   int sample_size = predictors.nrow();
   NumericVector norm_vector(sample_size);
   NumericVector knot(knot_size);
   double norm_xi, interval;
   int i, k;
   for (i = 0; i < sample_size; i++)
   {
      norm_xi = norm_2(predictors(i, _));
      norm_vector[i] = norm_xi;
   }
   double max_range = max(norm_vector);
   double increment = (2 * max_range) / double(knot_size + 1);
   for (k = 0; k < knot_size; k++)
   {
      interval = increment * double(k + 1);
      knot[k] = -max_range + interval;
   }
   knot.push_front(knot[0] - increment);
   knot.push_back(knot[knot_size] + increment);
   return knot;
}
//=========================================================================================
// norm_2
// [[Rcpp::export]]
double norm_2(NumericVector u)
{
   return sqrt(norm_square(u));
}
//=========================================================================================
// norm_square
// [[Rcpp::export]]
double norm_square(NumericVector u)
{
   return sum(pow(u, 2.0));
}
//=========================================================================================
// [[Rcpp::export]]
double zlambda(NumericVector a,
               double b,
               double c,
               NumericVector d,
               double lambda)
{
   NumericVector delta_d(d.size() + 1);
   IntegerVector order_a(d.size());
   NumericVector order_d(d.size());
   double q_lambda, q_lambda_zstar, zstar;
   //
   order_a = bubble_order(a);
   order_d = subsetNumVec(d, order_a);
   delta_d = delta(order_d);
   // initialize zstar by zero
   zstar = 0;
   // initialize q_lambda_zstar by q_lambda when z = 0
   q_lambda_zstar = q_lambda_point(zstar, a, b, c, d, lambda);
   // enumerate candidates of zstar
   NumericVector z(delta_d.size() + d.size());
   for (int j = 0; j < delta_d.size(); j++)
      z[j] = c - lambda * delta_d[j] / b;
   for (int k = delta_d.size(); k < z.size(); k++)
      z[k] = a[k - delta_d.size()];
   // find zstar
   for (int i = 0; i < z.size(); i++)
   {
      q_lambda = q_lambda_point(z[i], a, b, c, d, lambda);
      if (q_lambda < q_lambda_zstar)
      {
         zstar = z[i];
         q_lambda_zstar = q_lambda;
      }
   }
   return zstar;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector delta(NumericVector d)
{
   NumericVector del(d.size() + 1);
   double sum_d;
   sum_d = sum(d);
   int j;
   double a = 0;
   del[0] = -sum_d;
   del[1] = sum_d;
   for (j = 0; j < (d.size()-1); j++)
   {
      a += d[j];
      del[j+2] = 2 * a - sum_d;
   }
   return del;
}
//=========================================================================================
// computes q_lambda =
// (b / 2) * (z - c)^2 + lambda * ([d_1 * |z - a_1|] + ... + [d_k * |z - a_k|])
// [[Rcpp::export]]
double q_lambda_point(double z,
                      NumericVector a,
                      double b,
                      double c,
                      NumericVector d,
                      double lambda)
{
   int j;
   int d_size(d.size());
   double q_lambda = 0.5 * b * (z - c) * (z - c);
   for (j = 0; j < d_size; j++)
      q_lambda += lambda * d[j] * (std::abs(z - a[j]));
   return q_lambda;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector lambdas_all_alpha(int number_lambdas_alpha,
                                double lambda_alpha_max, double epsilon_lambda)
{
   // compute all lambdas
   // log l_k = log l_1 + (K - k) (log l_K - log l_1)/(K - 1), k = 1, ..., K
   // observe log l_1 = log l_K and log l_K = log l_1
   NumericVector lambdas_all(number_lambdas_alpha);
   double lambda_min = epsilon_lambda * lambda_alpha_max;
   //lambda_min = 1e-3;
   // lambda_max/lambda_min = 1 / epsilon_lambda = ratio_max_min
   double ratio_max_min = 1.0 / epsilon_lambda;
   double div, exponent;
   div = number_lambdas_alpha - 1;
   for (double lambda_index = 0; lambda_index < number_lambdas_alpha; lambda_index++)
   {
      exponent = lambda_index / div;
      lambdas_all[lambda_index] = lambda_min * pow(ratio_max_min, exponent);
   }
   lambdas_all[0] = 0.0;
   return lambdas_all;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector lambdas_all_beta(int number_lambdas,
                               double lambda_beta_max, double epsilon_lambda)
{
   // compute all lambdas
   // log l_k = log l_1 + (K - k) (log l_K - log l_1)/(K - 1), k = 1, ..., K
   // observe log l_1 = log l_K and log l_K = log l_1
   NumericVector lambdas_all(number_lambdas);
   double lambda_min = epsilon_lambda * lambda_beta_max;
   // lambda_max/lambda_min = 1 / epsilon_lambda = ratio_max_min
   double ratio_max_min = 1.0 / epsilon_lambda;
   double div, exponent;
   div = number_lambdas - 1;
   for (double lambda_index = 0; lambda_index < number_lambdas; lambda_index++)
   {
      exponent = lambda_index / div;
      lambdas_all[lambda_index] = lambda_min * pow(ratio_max_min, exponent);
   }
   lambdas_all[0] = 0.0;
   return lambdas_all;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector alpha2xi(NumericVector alpha)
{
   // input : alpha
   // output : xi
   int pm1 = alpha.size();
   int p = pm1 + 1;
   NumericVector xi(p);
   xi.fill(1);
   // xi[0]
   for(int i = 0; i < pm1; i++)
      xi[0] *= sinpi(alpha[i] / PI);
   // xi[1] ~ xi[p - 2]
   for (int j = 1; j < pm1; j++)
   {
      for (int i = 0; i < pm1 - j; i++)
         xi[j] *= sinpi(alpha[i] / PI);
      xi[j] *= cospi(alpha[pm1 - j] / PI);
   }
   // xi[pm1]
   xi[pm1] = cospi(alpha[0] / PI);
   return xi;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector xi2alpha(NumericVector xi)
{
   int p = xi.size();
   int pm1 = p - 1;
   double c;
   NumericVector alpha(pm1);
   for(int m = 0; m < pm1 - 1; m++)
   {
      c = xi[p - m - 1];
      for (int j = 0; j < m; j++)
         c /= (sin(alpha[j]) + 1e-10);
      alpha[m] = acos(c);
      //Rcout << "c = " << c << std::endl;
      //Rcout << "alpha[m] = " << alpha[m] << std::endl;
   }
   // alpha[pm1]
   alpha[pm1 - 1] = atan4(xi[0], xi[1]);
   return alpha;
}
//=========================================================================================
// [[Rcpp::export]]
NumericVector deriv_xi_j(NumericVector xi, NumericVector alpha, int j)
{
   // output : derivative vector of xi for alpha_j
   int p = xi.size();
   int pm1 = p - 1;
   int m, k;
   // index match
   j = j + 1;
   NumericVector deriv_xi(p);
   // 1
   deriv_xi[0] = 1;
   for(k = 0; k < j - 1; k++)
      deriv_xi[0] *= sin(alpha[k]);
   deriv_xi[0] *= cos(alpha[j - 1]);
   for(int k = j; k < pm1; k++)
      deriv_xi[0] *= sin(alpha[k]);
   // 2, ..., p - j
   for(m = 1; m < p - j; m++)
   {
      deriv_xi[m] = 1;
      for(k = 0; k < j - 1; k++)
         deriv_xi[m] *= sin(alpha[k]);
      deriv_xi[m] *= cos(alpha[j - 1]);
      for(int k = j; k < pm1 - m - 1; k++)
         deriv_xi[m] *= sin(alpha[k]);
      deriv_xi[m] *= cos(alpha[pm1 - m - 1]);
   }
   // p - j + 1
   deriv_xi[p - j] = 1;
   for(k = 0; k < j - 2; k++)
      deriv_xi[p - j] *= sin(alpha[k]);
   deriv_xi[p - j] *= -sin(alpha[j - 1]);
   return deriv_xi;
}
//=========================================================================================
// atan4 computes computes the inverse tangent of the ratio y/x.
// If y is not zero, then the tangent is y/x.
// Input
//    y, x = two quantities which represent the tangent of an angle.
// Output
//    an angle between 0 and 2 * pi, whose tangent is y/x,
//    and which lies in the appropriate quadrant so that the signs
//    of its cosine and sine match those of X and Y.
// [[Rcpp::export]]
double atan4(double y, double x)
{
   // local variables
   double abs_x, abs_y, theta;
   // output variable
   double atan_yx = 0.0;
   // Special cases:
   if (x == 0.0)
   {
      if (0.0 < y)
         atan_yx = PI / 2.0;
      else if (y < 0.0 )
         atan_yx = 3.0 * PI / 2.0;
      else if ( y == 0.0 )
         atan_yx = 0.0;
   }
   else if (y == 0.0)
   {
      if (0.0 < x)
         atan_yx = 0.0;
      else if (x < 0.0)
         atan_yx = PI;
   }
   //  We assume that atan2 is correct when both arguments are positive.
   else
   {
      abs_y = std::abs(y);
      abs_x = std::abs(x);
      theta = atan(abs_y / abs_x);
      if (0.0 < x && 0.0 < y)
         atan_yx = theta;
      else if (x < 0.0 && 0.0 < y)
         atan_yx = PI - theta;
      else if (x < 0.0 && y < 0.0)
         atan_yx = PI + theta;
      else if (0.0 < x && y < 0.0)
         atan_yx = 2.0 * PI - theta;
   }
   return atan_yx;
}
//=========================================================================================
// [[Rcpp::export]]
IntegerVector test_11(int a)
{
   return seq(0, a);
}
