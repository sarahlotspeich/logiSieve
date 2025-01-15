// [[Rcpp::depends(RcppEigen)]]

#ifndef UTILITY_H
#define UTILITY_H

#include <string>
#include <RcppEigen.h>

// Packages and function alias for speed testing
// #include <chrono>
// #include <ctime>
// const auto tic = std::chrono::system_clock::now;


using namespace std;
using namespace Eigen;
using namespace Rcpp;


/**
 * Prints an error message 'reason' and stops the program
 */
void stdError (const string reason);

/**
 * Returns true if arr1 has a larger value than arr2 closer to index 0
 * Returns false if equal or arr2 has the larger value
 */
bool BigArray (const RowVectorXd& arr1, const RowVectorXd& arr2, const int& maxIndex);

/**
 * Returns true if arr1 has a smaller value than arr2 closer to index 0
 * Returns false if equal or arr2 has the smaller value
 */
bool SmallArray (const RowVectorXd& arr1, const RowVectorXd& arr2, const int& maxIndex);

/** 
 *  Check if two arrays are equal to a certain index
 */
bool EqualArray (const RowVectorXd& arr1, const RowVectorXd& arr2, const int& maxIndex);

/**
 *
 */
VectorXi indexx_Matrix_Row (const MatrixXd& mat);

/**
 *  Find the number of unique rows in "mat" and return it as "n1". "index" is the index
 *  vector of the rows of "mat" output by "indexx_Matrix_Row". "mat" and "index" are not
 *  changed.
 */
int Num_Uni_Matrix_Row (const MatrixXd& mat, const VectorXi& index) ;

/**
 *  Create unique covariates. "mat" is a matrix (not changed). "index" is the index vector
 *  of the rows of mat output by indexx_Matrix_Row (not changed). "uni_mat" is the matrix
 *  of unique observations (output). "ind" is the vector of indexes of which one of the
 *  unique rows the original observations correspond to (output). "N_uni" stores the
 *  number of observations for each distinct covariates row (output).
 */
void Create_Uni_Matrix_Row (const MatrixXd& mat, const VectorXi& index, MatrixXd& uni_mat, VectorXi& ind, VectorXi& N_uni) ;

/**
 * 
 */
VectorXi indexx_Vector (const VectorXd& vec) ;

/***************************************************************************************
 Find the number of unique events, and return it as "n_event".
 "Y" is the vector of times.
 "Y_index" is the index vector of Y output by "indexx_Vector()".
 "Delta" is the vector of event indicators.
***************************************************************************************/
int Num_Distinct_Events (const VectorXd& Y, const VectorXi& Y_index, const VectorXi& Delta) ;

/***************************************************************************************
 Create unique events.
 "Y" is the vector of times.
 "Y_index" is the index vector of Y output by "indexx_Vector()".
 "Delta" is the vector of event indicators.
 "Y_uni_event" is a vector of unique events.
 "Y_risk_ind" is a vector of indexes to which one of the risk sets each element in "Y" corresponds.
 "Y_uni_event_n" is a vector of numbers of events corresponding to each unique event time.
***************************************************************************************/
void Create_Uni_Events (const VectorXd& Y, const VectorXi& Y_index, const VectorXi& Delta, 
    VectorXd& Y_uni_event, VectorXi& Y_risk_ind, VectorXi& Y_uni_event_n); 

/***************************************************************************************
 Create unique events.
 "Y" is the vector of times.
 "Y_index" is the index vector of Y output by "indexx_Vector()".
 "Delta" is the vector of event indicators.
 "Y_uni_event" is a vector of unique events.
 "Y_risk_ind" is a vector of indexes to which one of the risk sets each element in "Y" corresponds.
 "Y_uni_event_n" is a vector of numbers of events corresponding to each unique event time.
 "L" is the vector of left-truncation times.
 "L_index" is the index vector of L output by "indexx_Vector()".
 "L_risk_ind" is a vector of indexes to which one of the risk sets each element in "L" corresponds.
***************************************************************************************/
void Create_Uni_Events_LeftTrunc (const VectorXd& Y, const VectorXd& L, const VectorXi& Y_index, const VectorXi& L_index,
    const VectorXi& Delta, VectorXd& Y_uni_event, VectorXi& Y_risk_ind, VectorXi& Y_uni_event_n, VectorXi& L_risk_ind) ;

/**
 * Calculate the variance of an array via difference of squares (R^2)
 */
double Var(const VectorXd& arr);

#endif
