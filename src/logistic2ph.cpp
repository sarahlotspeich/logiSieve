// This definition allows us to do some big matrix multiplication (calculateHessian)
// Set in Makevars
// #define ARMA_64BIT_WORD 1

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <stdlib.h>


using namespace Rcpp;
using namespace std;

//' Multiply matrix by vector
//'
//' Multiplies each column of a matrix by a vector
//'
//' @param mat The matrix
//' @param v The vector
//' @return mat * v
//' @noRd
arma::mat matTimesVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  // Never an issue in this project, but if you're copy+pasting this code into yours,
  // you may want to uncomment the error checking
  // if (v.n_elem < mat.n_rows)
  // {
  //   arma::vec oldV = v;
  //   for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
  //   {
  //     v = join_vert(v, oldV);
  //   }
  // }

  // Multiply each col by the vector
  mat.each_col() %= v;

  return mat;
}

//' Divide matrix by vector
//'
//' Divides each column of a matrix by a vector
//'
//' @param mat The matrix
//' @param v The vector divisor
//' @return mat / v, or mat * v^-1
//' @noRd
arma::mat matDivideVec(arma::mat mat, arma::vec v)
{
  // Ensure the vector is the right length
  if (v.n_elem < mat.n_rows)
  {
    arma::vec oldV = v;
    for (int i = 1; i < (int)mat.n_rows / (int)oldV.n_elem; ++i )
    {
      v = join_vert(v, oldV);
    }
  }

  // Divide each col by the vector
  mat.each_col() /= v;

  return mat;
}


// TRANSLATING PACKAGE FUNCTIONS TO CPP FOR SPEED BOOST

//' Prepend ones to a w_t
//'
//' Lengthens a vector by prepending n ones
//'
//' @param w_t_original The original vector
//' @param n The number of ones to add to the front of the vector
//' @param modifyW_T If false, instantly returns w_t_original
//' @noRd
// [[Rcpp::export(.lengthenWT)]]
arma::vec lengthenWT(
  const arma::vec& w_t_original,
  const int& n,
  const bool& modifyW_T = true)
{
  // If we don't modify w_t_original, return it immediately
  if (!modifyW_T)
    return w_t_original;

  // If w_t_original is empty, return the ones vector of length n
  if (w_t_original.is_empty())
      return arma::ones<arma::vec>(n);

  // Put n 1's in front of w_t_original
  arma::vec w_t(n + w_t_original.n_elem);   // Initialize vector
  w_t = join_vert(arma::ones<arma::vec>(n), w_t_original);  // Join vector of 1's with w_t_original


  return w_t;
}

//' Calculate Mu
//'
//' Calculates the value of mu according to two variables
//' Small helper function
//'
//' @param design_mat The design matrix
//' @param prev The previous iteration of the design matrix
//' @noRd
// [[Rcpp::export(.calculateMu)]]
arma::vec calculateMu(
  const arma::mat& design_mat,
  const arma::mat& prev)
{
  arma::mat mu = (design_mat * prev) * -1;
  arma::vec mu1 = exp(mu).as_col();
  return mu1 / (1 + mu1);
}

//' Calculate gradient
//'
//' Calculates a gradient given w_t and a design matrix
//' TODO
//'
//' @param w_t A vector indicating ??
//' @param n The number of ones to prepend to w_t
//' @param design_mat The design matrix
//' @param Y_col The column of validated Y values from the complete data matrix
//' @param muVector The vector calculated by calculateMu
//' @param modifyW_T Whether to add ones to the beginning of w_t
//' @noRd
// [[Rcpp::export(.calculateGradient)]]
arma::vec calculateGradient(
  arma::vec& w_t,
  const int& n,
  const arma::mat& design_mat,
  const arma::vec& Y_col,
  const arma::vec& muVector,
  const bool& modifyW_T = false)
{
  // Put n 1's in front of w_t
  w_t = lengthenWT(w_t, n, modifyW_T);

  // Calculate gradient
  arma::vec sumsVector = Y_col - 1 + muVector;

  // Convert to NumericMatrix for Rcpp sugar's colSums
  arma::mat temp = matTimesVec(design_mat, sumsVector);
  temp = matTimesVec(temp, w_t);

  // Sum's default behavior on matrices is like colSums
  arma::rowvec gradient = sum(temp);


  return reshape(gradient, gradient.n_elem, 1);

}

//' Calculate Hessian Matrix
//'
//' Calculates the Hessian Matrix and lengthens w_t by n
//'
//' @param design_mat The design matrix
//' @param w_t The vector ??
//' @param muVector The vector returned by calculateMu
//' @param n The number of ones to prepend to w_t
//' @param mus An empty, pre-allocated vector of the same length as muVector, pre-allocated memory saves time
//' @param modifyW_T Whether to add ones to the beginning of w_t
//' @noRd
// [[Rcpp::export(.calculateHessian)]]
arma::mat calculateHessian(
  const arma::mat& design_mat,
  arma::vec& w_t,
  const arma::vec& muVector,
  const int & n,
  arma::vec& mus,
  const bool& modifyW_T = false  )
{
  w_t = lengthenWT(w_t, n, modifyW_T);

  // post_multiply = c(w_t * muVector * (muVector - 1)) * gamma_design_mat
  mus = muVector % (muVector - 1);
  mus = mus % w_t;

  return design_mat.t() * matTimesVec(design_mat, mus);

}

//' Calculate pYstar
//'
//' TODO
//'
//' @param gamma_design_mat The gamma design matrix
//' @param startRow The number of the rows to start within gamma_design_mat
//' @param prev_gamma The previous iteration of gamma_design_mat
//' @param comp_dat_all The complete dataset
//' @param Y_unval_index Which column of comp_dat_all houses the unvalidated Y variable
//' @noRd
// [[Rcpp::export(.pYstarCalc)]]
arma::vec pYstarCalc(
  const arma::mat& gamma_design_mat,
  const int& startRow,
  const arma::mat& prev_gamma,
  const arma::mat& comp_dat_all,
  const int& Y_unval_index)
{
  // pYstar and mu_gamma are pre-allocated to save time with memory management

  // same as gamma_design_mat[-c(1:n),]
  // get the elements of gamma_design_mat excluding the first excludeRows rows
  arma::mat filtered_gamma_design_mat = gamma_design_mat.rows(startRow-1, gamma_design_mat.n_rows-1);
  arma::vec mu_gamma = filtered_gamma_design_mat * prev_gamma;
  arma::vec pYstar = 1 / (1 + exp(mu_gamma * -1));

  arma::vec checkVector = comp_dat_all.col(Y_unval_index).rows(startRow-1, comp_dat_all.n_rows-1);
  for (unsigned int i = 0; i < pYstar.size(); ++i)
  {
    if (checkVector(i) == 0)
    {
      pYstar(i) = 1 - pYstar(i);
    }
  }

  return pYstar;
}

