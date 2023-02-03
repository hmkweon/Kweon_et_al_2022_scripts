#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.hpp"
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

Eigen::MatrixXd Standardize(const Eigen::MatrixXd &mat)
{
    Eigen::MatrixXd std_mat(mat.rows(), mat.cols());
    for (int i = 0; i < mat.cols(); ++i)
    {
        Eigen::ArrayXd r = mat.col(i).array();
        double colMean = r.mean();
        double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
        std_mat.col(i) = (r - colMean) / colSdev;
    }
    return std_mat;
}

// [[Rcpp::export]]
Eigen::MatrixXd Stdz(const Eigen::Map<Eigen::MatrixXd> &mat)
{
    Eigen::MatrixXd std_mat(mat.rows(), mat.cols());
    for (int i = 0; i < mat.cols(); ++i)
    {
        Eigen::ArrayXd r = mat.col(i).array();
        double colMean = r.mean();
        double colSdev = sqrt((r - colMean).square().sum() / (mat.rows() - 1));
        std_mat.col(i) = (r - colMean) / colSdev;
    }
    return std_mat;
}


int count_if(LogicalVector x)
{
    int counter = 0;
    for (int i = 0; i < x.size(); i++)
    {
        if (x[i] == TRUE)
        {
            counter++;
        }
    }
    return counter;
}

Eigen::MatrixXd submat(const Eigen::MatrixXd &X, const LogicalVector &mask)
{ // this function converts data frame to matrix while subsetting
    int nRows = count_if(mask);
    NumericMatrix mat(nRows, X.size());
    for (int i = 0; i < X.cols(); i++)
    {
        NumericVector new_col = wrap(X.col(i));
        new_col = new_col[mask];
        mat(_, i) = new_col;
    }
    return as<Eigen::Map<Eigen::MatrixXd>>(mat);
}


Eigen::MatrixXd DFtoNM(const DataFrame &X, const IntegerVector &IND)
{ // this function converts data frame to matrix while subsetting
    int nRows = IND.size();
    NumericMatrix mat(nRows, X.size());
    for (int i = 0; i < X.size(); i++)
    {
        NumericVector new_col = X[i];
        new_col = new_col[IND - 1];
        mat(_, i) = new_col;
    }
    return as<Eigen::Map<Eigen::MatrixXd>>(mat);
}


Eigen::MatrixXd DFtoNM2(const DataFrame &X)
{ // this function converts data frame to matrix
    NumericMatrix mat(X.nrows(), X.size());
    for (int i = 0; i < X.size(); i++)
    {
        NumericVector new_col = X[i];
        mat(_, i) = new_col;
    }
    return as<Eigen::Map<Eigen::MatrixXd>>(mat);
}

// [[Rcpp::export]]
Eigen::MatrixXd make_res(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y)
{ // resdiaulize Y

    Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(X.cols(), X.cols()).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint()));
    Eigen::MatrixXd E = (Y - X * llt.solve(X.transpose() * Y));

    return E;
}

// These functions will be slower if n >> p

// [[Rcpp::export]]
Eigen::MatrixXd get_XXt(const Eigen::Map<Eigen::MatrixXd> &X)
{
    return Eigen::MatrixXd(X.rows(), X.rows()).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X);
}

Eigen::MatrixXd AAt(const Eigen::MatrixXd &X)
{
    return Eigen::MatrixXd(X.rows(), X.rows()).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X);
}

// using selfadjointview slower if p>>n i.e., numer of vars >>> number of indiv.

// [[Rcpp::export]]
Eigen::MatrixXd get_XtX(const Eigen::Map<Eigen::MatrixXd> &X)
{
    return Eigen::MatrixXd(X.cols(), X.cols()).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint());
}

Eigen::MatrixXd AtA(const Eigen::MatrixXd &X)
{
    return Eigen::MatrixXd(X.cols(), X.cols()).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint());
}

// [[Rcpp::export]]
Eigen::MatrixXd get_cholinv(const Eigen::Map<Eigen::MatrixXd> &X)
{
    int p = X.cols();
    Eigen::LLT<Eigen::MatrixXd> llt_denom((X).selfadjointView<Eigen::Upper>());

    return llt_denom.solve(Eigen::MatrixXd::Identity(p, p));
}

\
