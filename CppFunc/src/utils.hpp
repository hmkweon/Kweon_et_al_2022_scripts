#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

Eigen::MatrixXd Standardize(const Eigen::MatrixXd &mat);

int count_if(LogicalVector x);

Eigen::MatrixXd submat(const Eigen::MatrixXd &X, const LogicalVector &mask);

Eigen::MatrixXd make_res(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y);

Eigen::MatrixXd DFtoNM(const DataFrame &X, const IntegerVector &IND);

Eigen::MatrixXd DFtoNM2(const DataFrame &X);

Eigen::MatrixXd AtA(const Eigen::MatrixXd &X);

Eigen::MatrixXd AAt(const Eigen::MatrixXd &X);

List householderQR(const Eigen::Map<Eigen::MatrixXd> &X);
#endif