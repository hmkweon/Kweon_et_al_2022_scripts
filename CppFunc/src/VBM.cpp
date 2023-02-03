#include <Rcpp.h>
#include <RcppEigen.h>
#include "utils.hpp"
//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;



// [[Rcpp::export]]
Eigen::MatrixXd reg_F(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const DataFrame &Y, const IntegerVector &IND)
{

    Eigen::MatrixXd stdY = Standardize(DFtoNM(Y, IND));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(3).bottomRows(2); //return only beta for SES variables
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(B.rows() + 2, B.cols());
    output << B, SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd reg_F2(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const DataFrame &Y)
{
    // the version that doesn't need subset individuals
    Eigen::MatrixXd stdY = Standardize(DFtoNM2(Y));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(3).bottomRows(2); //return only beta for SES variables
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(B.rows() + 2, B.cols());
    output << B, SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd reg_F3(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const DataFrame &Y, const int &K)
{
    // the version that doesn't need subset individuals
    Eigen::MatrixXd stdY = Standardize(DFtoNM2(Y));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(K + 1).bottomRows(K); //return first K coeffs
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(B.rows() + 2, B.cols());
    output << B, SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd reg_F4(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const DataFrame &Y, const IntegerVector &IND, const int &K)
{

    Eigen::MatrixXd stdY = Standardize(DFtoNM(Y, IND));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(K + 1).bottomRows(K); //return first K coeffs
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(B.rows() + 2, B.cols());
    output << B, SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd reg_cov(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, const DataFrame &Y, const IntegerVector &IND, const int &K)
{
    // returns sigma
    int p = R.cols();

    Eigen::MatrixXd stdY = Standardize(DFtoNM(Y, IND));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(K + 1).bottomRows(K); //return first K coeffs

    Eigen::MatrixXd sigma_sq = (stdY - Q * QtY).colwise().squaredNorm() / (Y.rows() - p);

    Eigen::MatrixXd output(B.rows() + 1, B.cols());
    output << B, sigma_sq;
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd reg_F_ROI(const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, const Eigen::Map<Eigen::MatrixXd> &Q_r, const Eigen::Map<Eigen::MatrixXd> Y, const int &K)
{

    Eigen::MatrixXd stdY = Standardize(Y);
    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(K + 1).bottomRows( K ); //return only beta for SES variables
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(B.rows() + 2, B.cols());
    output << B, SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd res_vox(const Eigen::Map<Eigen::MatrixXd> &Q, const DataFrame &Y, const IntegerVector &IND)
{

    Eigen::MatrixXd stdY = Standardize(DFtoNM(Y, IND));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;

    return (stdY - Q * QtY);
}

// [[Rcpp::export]]
Eigen::MatrixXd res_ROI(const Eigen::Map<Eigen::MatrixXd> &Q, const DataFrame &Y)
{

    Eigen::MatrixXd stdY = Standardize(DFtoNM2(Y));
    Eigen::MatrixXd QtY = Q.transpose() * stdY;

    return (stdY - Q * QtY);
}

// [[Rcpp::export]]
Eigen::MatrixXd perm_mri(const Eigen::Map<Eigen::MatrixXd> &Q, const Eigen::Map<Eigen::MatrixXd> &Q_r, const Eigen::Map<Eigen::MatrixXd> Y, const int DF, const int K)
{

    Eigen::ArrayXd SSR = (Y - Q * (Q.transpose() * Y)).colwise().squaredNorm();
    Eigen::ArrayXd SSR_r = (Y - Q_r * (Q_r.transpose() * Y)).colwise().squaredNorm();

    return (((SSR_r - SSR) / K) / (SSR / DF)).matrix();
}


// [[Rcpp::export]]
Eigen::MatrixXd perm_mri_t(const Eigen::Map<Eigen::MatrixXd> &Q, const Eigen::Map<Eigen::MatrixXd> &Q_r, const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Y, const Eigen::Map<Eigen::VectorXd> &varX_div_DF)
{
    int K = varX_div_DF.size();
    Eigen::MatrixXd QtY = Q.transpose() * Y;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY).topRows(K + 1).bottomRows(K); //return only beta for SES variables

    Eigen::MatrixXd SSR = (Y - Q * QtY).colwise().squaredNorm();
    Eigen::MatrixXd SE = (varX_div_DF * SSR).cwiseSqrt();

    Eigen::MatrixXd SSR_r = (Y - Q_r * (Q_r.transpose() * Y)).colwise().squaredNorm();
    Eigen::MatrixXd F = ((SSR_r - SSR) / K).cwiseProduct( (SSR/ (Q.rows() - R.rows())).cwiseInverse() );

    Eigen::MatrixXd OUT(B.rows() + F.rows(), B.cols());
    OUT << B.cwiseProduct(SE.cwiseInverse()), F;
    
    return OUT;
}

// [[Rcpp::export]]
Eigen::MatrixXd IV_F(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const DataFrame &Y, const IntegerVector &IND, const int &K)
{

    Eigen::MatrixXd stdY = Standardize(DFtoNM(Y, IND));

    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY);

    Eigen::MatrixXd SSR_2sls = (stdY - X * B).colwise().squaredNorm();
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(K + 3, B.cols());
    output << B.topRows(K + 1).bottomRows(K), SSR_2sls.transpose(), SSR.transpose(), SSR_r.transpose();
    return output;
}

// [[Rcpp::export]]
Eigen::MatrixXd IV_F_ROI(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &R, const Eigen::Map<Eigen::MatrixXd> &Q, Eigen::Map<Eigen::MatrixXd> &Q_r, const Eigen::Map<Eigen::MatrixXd> &Y, const int &K)
{

    Eigen::MatrixXd stdY = Standardize(Y);

    Eigen::MatrixXd QtY = Q.transpose() * stdY;
    Eigen::MatrixXd B = R.triangularView<Eigen::Upper>().solve(QtY);

    Eigen::MatrixXd SSR_2sls = (stdY - X * B).colwise().squaredNorm();
    Eigen::VectorXd SSR = (stdY - Q * QtY).colwise().squaredNorm();
    Eigen::VectorXd SSR_r = (stdY - Q_r * (Q_r.transpose() * stdY)).colwise().squaredNorm();

    Eigen::MatrixXd output(K + 3, B.cols());
    output << B.topRows(K + 1).bottomRows(K), SSR_2sls.transpose(), SSR.transpose(), SSR_r.transpose();
    return output;
}