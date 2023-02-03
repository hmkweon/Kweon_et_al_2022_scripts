#include <Rcpp.h>
#include <RcppEigen.h>
#include <Rcpp/Benchmark/Timer.h>
#include "utils.hpp"

//[[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
List lm_cpp(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y)
{

  int p = X.cols();
  Eigen::MatrixXd XtY = X.transpose() * Y;
  Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint()));
  Eigen::MatrixXd B = llt.solve(XtY);

  // double sigma = (Y - X * B).squaredNorm() / (Y.rows() - p);
  // Eigen::VectorXd SE = sqrt(sigma) * (llt.solve(Eigen::MatrixXd::Identity(p, p))).diagonal().cwiseSqrt();
  Eigen::MatrixXd sigma = ((Y - X * B).colwise().squaredNorm() / (Y.rows() - p)).cwiseSqrt();
  Eigen::MatrixXd SE = (llt.solve(Eigen::MatrixXd::Identity(p, p))).diagonal().cwiseSqrt() * sigma;

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
  // Named("R2") = 1 - sigma * (Y.rows() - p) / (Y.array() - Y.array().mean()).matrix().squaredNorm() );
}

// [[Rcpp::export]]
List lm_cpp_het(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y)
{

  int p = X.cols();
  Eigen::MatrixXd XtY = X.transpose() * Y;
  Eigen::LLT<Eigen::MatrixXd> llt_XtX(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint()));
  Eigen::MatrixXd B = llt_XtX.solve(XtY);

  Eigen::MatrixXd INV = llt_XtX.solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd U = (Y - X * B);
  Eigen::MatrixXd SE(X.cols(), Y.cols());
  for (int i = 0; i < SE.cols(); ++i)
  {
    SE.col(i) = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(INV.selfadjointView<Eigen::Upper>() * (X.transpose() * U.col(i).cwiseAbs().asDiagonal())).diagonal().cwiseSqrt();
  }

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
  // Named("R2") = 1 - sigma * (Y.rows() - p) / (Y.array() - Y.array().mean()).matrix().squaredNorm() );
}

// [[Rcpp::export]]
List lm_cpp_cluster(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y, const NumericVector &C)
{

  int p = X.cols();
  Eigen::MatrixXd XtY = X.transpose() * Y;
  Eigen::LLT<Eigen::MatrixXd> llt_XtX(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(X.adjoint()));
  Eigen::MatrixXd B = llt_XtX.solve(XtY);

  NumericVector C_u = unique(C);

  Eigen::MatrixXd INV = llt_XtX.solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd U = (Y - X * B);
  Eigen::MatrixXd SE(X.cols(), Y.cols());
  Eigen::MatrixXd Sigma_temp;
  Eigen::MatrixXd VCOV_temp;

  for (int i = 0; i < SE.cols(); ++i)
  {
    NumericVector U_i = wrap(U.col(i));
    Eigen::MatrixXd XtE(Eigen::MatrixXd(p, p).setZero());

    for (int c = 0; c < C_u.size(); ++c)
    {
      // make mask
      NumericVector c_vec(C.size());
      c_vec.fill(C_u[c]);
      LogicalVector mask = (C != c_vec);
        // try change this, using .block
      XtE.selfadjointView<Eigen::Upper>().rankUpdate(submat(X, mask).transpose() * as<Eigen::Map<Eigen::VectorXd>>(U_i[mask]));
    }

    SE.col(i) = (INV.selfadjointView<Eigen::Upper>() * (XtE.selfadjointView<Eigen::Upper>() * INV)).diagonal().cwiseSqrt();
  }

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
 
}


// [[Rcpp::export]]
List iv_cpp(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Z, const Eigen::Map<Eigen::MatrixXd> &Y)
{

  int p = X.cols();

  Eigen::LLT<Eigen::MatrixXd> llt_ZtZ(AtA(Z));
  // Eigen::MatrixXd PZ = Z * llt.solve(Eigen::MatrixXd::Identity(p, p)).selfadjointView<Eigen::Upper>() * Z.transpose();
  Eigen::MatrixXd Uinv = llt_ZtZ.matrixU().solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd Z_Uinv = Z * Uinv.triangularView<Eigen::Upper>();
  Eigen::MatrixXd Xt_Z_Uinv = X.transpose() * Z_Uinv;
  Eigen::MatrixXd Xt_PZ_X = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Xt_Z_Uinv);

  Eigen::LLT<Eigen::MatrixXd> llt_denom((Xt_PZ_X).selfadjointView<Eigen::Upper>());
  Eigen::MatrixXd B = llt_denom.solve(Xt_Z_Uinv * (Z_Uinv.transpose() * Y));
  Eigen::MatrixXd sigma = ((Y - X * B).colwise().squaredNorm() / (Y.rows() - p)).cwiseSqrt();
  Eigen::MatrixXd SE = llt_denom.solve(Eigen::MatrixXd::Identity(p, p)).diagonal().cwiseSqrt() * sigma;

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
}

// [[Rcpp::export]]
List iv_cpp_het(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Z, const Eigen::Map<Eigen::MatrixXd> &Y)
{

  int p = X.cols();

  Eigen::LLT<Eigen::MatrixXd> llt_ZtZ(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Z.adjoint()));
  // Eigen::MatrixXd PZ = Z * llt.solve(Eigen::MatrixXd::Identity(p, p)).selfadjointView<Eigen::Upper>() * Z.transpose();
  Eigen::MatrixXd Uinv = llt_ZtZ.matrixU().solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd Z_Uinv = Z * Uinv.triangularView<Eigen::Upper>();
  Eigen::MatrixXd Xt_Z_Uinv = X.transpose() * Z_Uinv;
  Eigen::MatrixXd Xt_PZ_X = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Xt_Z_Uinv);

  Eigen::LLT<Eigen::MatrixXd> llt_denom((Xt_PZ_X).selfadjointView<Eigen::Upper>());
  Eigen::MatrixXd B = llt_denom.solve(Xt_Z_Uinv * (Z_Uinv.transpose() * Y));
  Eigen::MatrixXd U = (Y - X * B);
  Eigen::MatrixXd INV = llt_denom.solve(Eigen::MatrixXd::Identity(p, p));

  Eigen::MatrixXd SE(X.cols(), Y.cols());
  for (int i = 0; i < SE.cols(); ++i)
  {
    SE.col(i) = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(INV.selfadjointView<Eigen::Upper>() * Xt_Z_Uinv * (Z_Uinv.transpose() * U.col(i).cwiseAbs().asDiagonal())).diagonal().cwiseSqrt();
  }

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
}

// [[Rcpp::export]]
List iv_cpp_cluster(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Z, const Eigen::Map<Eigen::MatrixXd> &Y, const NumericVector &C)
{

  int p = X.cols();

  Eigen::LLT<Eigen::MatrixXd> llt_ZtZ(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Z.adjoint()));
  // Eigen::MatrixXd PZ = Z * llt.solve(Eigen::MatrixXd::Identity(p, p)).selfadjointView<Eigen::Upper>() * Z.transpose();
  Eigen::MatrixXd Uinv = llt_ZtZ.matrixU().solve(Eigen::MatrixXd::Identity(p, p));
  Eigen::MatrixXd Z_Uinv = Z * Uinv.triangularView<Eigen::Upper>();
  Eigen::MatrixXd Xt_Z_Uinv = X.transpose() * Z_Uinv;
  Eigen::MatrixXd Xt_PZ_X = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Xt_Z_Uinv);

  Eigen::LLT<Eigen::MatrixXd> llt_denom((Xt_PZ_X).selfadjointView<Eigen::Upper>());
  Eigen::MatrixXd B = llt_denom.solve(Xt_Z_Uinv * (Z_Uinv.transpose() * Y));
  Eigen::MatrixXd U = (Y - X * B);
  Eigen::MatrixXd INV = llt_denom.solve(Eigen::MatrixXd::Identity(p, p));
  NumericVector C_u = unique(C);

  Eigen::MatrixXd SE(X.cols(), Y.cols());
  for (int i = 0; i < SE.cols(); ++i)
  {
    NumericVector U_i = wrap(U.col(i));
    Eigen::MatrixXd Sigma_i = Eigen::MatrixXd(U_i.size(), U_i.size()).setZero().selfadjointView<Eigen::Upper>();

    for (int c = 0; c < C_u.size(); ++c)
    {
      NumericVector U_i_c = clone(U_i);

      //subset U_i for cluster c
      NumericVector c_vec(C.size());
      c_vec.fill(C_u[c]);
      LogicalVector mask = (C != c_vec);
      U_i_c[mask] = 0;

      Sigma_i.selfadjointView<Eigen::Upper>().rankUpdate(as<Eigen::Map<Eigen::VectorXd>>(U_i_c));
    }

    SE.col(i) = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(INV.selfadjointView<Eigen::Upper>() * Xt_Z_Uinv * (Z_Uinv.transpose() * Sigma_i.llt().matrixL())).diagonal().cwiseSqrt();
  }

  return List::create(Named("BETA") = B,
                      Named("SE") = SE);
}


// [[Rcpp::export]]
Eigen::MatrixXd RIDGE_K(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y, const int &K)
{
  // Eigen::MatrixXd Y = Standardize(Y_input);
  // Eigen::MatrixXd X = Standardize(X_input);

  int p = X.cols();
  Eigen::MatrixXd XtX = AtA(X);
  
  Eigen::LLT<Eigen::MatrixXd> denom_llt((XtX + (Eigen::MatrixXd::Identity(p, p) * K)).selfadjointView<Eigen::Upper>());
  Eigen::MatrixXd B = denom_llt.solve(X.transpose() * Y);

  return B;
}


Eigen::ArrayXd Dplus(const Eigen::ArrayXd &d)
{
  Eigen::ArrayXd di(d.size());
  for (int j = 0; j < d.size(); ++j)  di[j] = 1 / d[j];
  return di;
}

// [[Rcpp::export]]
Eigen::MatrixXd RIDGE_multi_K(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y, const Eigen::Map<Eigen::VectorXd> &K)
{

      Eigen::BDCSVD<Eigen::MatrixXd>
          UDV(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

  int p = X.cols();

  // this needs to be computed once.
  Eigen::MatrixXd DUtY = UDV.singularValues().matrix().asDiagonal() * (UDV.matrixU().transpose() * Y);
  // Eigen::MatrixXd VDsq_p = UDV.matrixV() * (Dplus(UDV.singularValues())).pow(2).matrix().asDiagonal();

  Eigen::MatrixXd B(p, K.size());  
  for (int i = 0; i < K.size(); ++i)
    {
      Eigen::MatrixXd VDsq_p = UDV.matrixV() * (Dplus(UDV.singularValues().array().pow(2) + K[i])).matrix().asDiagonal();
      B.col(i) = VDsq_p.matrix() * DUtY;
      // B.col(i) = (VDsq_p.array()/K[i]).matrix() * DUtY;
    }

    return B;
}

// [[Rcpp::export]]
Eigen::MatrixXd RIDGE_L(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y, const Eigen::Map<Eigen::VectorXd> &L)
{
  // Eigen::MatrixXd Y = Standardize(Y_input);
  // Eigen::MatrixXd X = Standardize(X_input);

  int p = X.cols();
  Eigen::MatrixXd XtX = AtA(X);
  Eigen::MatrixXd XtY = X.transpose() * Y;

  Eigen::MatrixXd B(p, L.size());
  for (int i = 0; i < L.size(); ++i)
  {
    Eigen::LLT<Eigen::MatrixXd> denom_llt((XtX + (Eigen::MatrixXd::Identity(p, p) * L[i])).selfadjointView<Eigen::Upper>());
    B.col(i) = denom_llt.solve(XtY);
  }
  return B;
}


// [[Rcpp::export]]
List RIDGE_L_AIC(const Eigen::Map<Eigen::MatrixXd> &X, const Eigen::Map<Eigen::MatrixXd> &Y, const Eigen::Map<Eigen::VectorXd> &L)
{
  // Eigen::MatrixXd Y = Standardize(Y_input);
  // Eigen::MatrixXd X = Standardize(X_input);

  int p = X.cols();
  Eigen::MatrixXd XtX = AtA(X);
  Eigen::MatrixXd XtY = X.transpose() * Y;

  Eigen::MatrixXd B(p, L.size());
  Eigen::VectorXd AIC(L.size());  
  Eigen::VectorXd DF(L.size());  
  double pi = 3.141592653589793238463;

  for (int i = 0; i < L.size(); ++i)
  {
    Eigen::LLT<Eigen::MatrixXd> denom_llt((XtX + (Eigen::MatrixXd::Identity(p, p) * L[i])).selfadjointView<Eigen::Upper>());
    B.col(i) = denom_llt.solve(XtY);

    double SSR = (Y - X * B).squaredNorm();
    DF(i) = ((X * denom_llt.solve(Eigen::MatrixXd::Identity(p, p)).selfadjointView<Eigen::Upper>() ) * X.transpose()).diagonal().sum();
    AIC(i) = 3*DF(i) + 2*Y.rows() * log( sqrt(2*pi*SSR / DF(i)))  ;
  }

  return List::create(Named("B") = B,
                      Named("AIC") = AIC,
                      Named("DF") = DF);
}


// [[Rcpp::export]]
Eigen::MatrixXd RIDGE_L_Kfolds(const Eigen::Map<Eigen::MatrixXd> &XtX, const Eigen::Map<Eigen::MatrixXd> &XtY, const Eigen::Map<Eigen::VectorXd> &L)
{

  int p = XtX.cols();

  Eigen::MatrixXd B(p, L.size());
  for (int i = 0; i < L.size(); ++i)
  {
    Eigen::LLT<Eigen::MatrixXd> denom_llt((XtX + (Eigen::MatrixXd::Identity(p, p) * L[i])).selfadjointView<Eigen::Upper>());
    B.col(i) = denom_llt.solve(XtY);
  }
  return B;
}

// [[Rcpp::export]]
List SVD_cpp(const Eigen::Map<Eigen::MatrixXd> &X)
{
  Eigen::BDCSVD<Eigen::MatrixXd> UDV(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

  return List::create(Named("d") = UDV.singularValues(),
                      Named("U") = UDV.matrixU(),
                      Named("V") = UDV.matrixV());
}

  // [[Rcpp::export]]
  List RIDGE_IV_K(const Eigen::Map<Eigen::MatrixXd> &X_input, const Eigen::Map<Eigen::MatrixXd> &Z_input, const Eigen::Map<Eigen::MatrixXd> &Y_input, const int &K)
  {

    Eigen::MatrixXd Y = Standardize(Y_input);
    Eigen::MatrixXd X = Standardize(X_input);
    Eigen::MatrixXd Z = Standardize(Z_input);
    int p = X.cols();

    Eigen::LLT<Eigen::MatrixXd> llt(Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Z.adjoint()));
    // Eigen::MatrixXd PZ = Z * llt.solve(Eigen::MatrixXd::Identity(p, p)).selfadjointView<Eigen::Upper>() * Z.transpose();
    Eigen::MatrixXd Uinv = llt.matrixU().solve(Eigen::MatrixXd::Identity(p, p));
    Eigen::MatrixXd Z_Univ = Z * Uinv.triangularView<Eigen::Upper>();
    Eigen::MatrixXd Xt_Z_Univ = X.transpose() * Z_Univ;
    Eigen::MatrixXd Xt_PZ_X = Eigen::MatrixXd(p, p).setZero().selfadjointView<Eigen::Upper>().rankUpdate(Xt_Z_Univ);

    Eigen::LLT<Eigen::MatrixXd> llt_denom((Xt_PZ_X + Eigen::MatrixXd::Identity(p, p) * K).selfadjointView<Eigen::Upper>());
    Eigen::VectorXd B = llt_denom.solve(Xt_Z_Univ * (Z_Univ.transpose() * Y));
    Eigen::MatrixXd INV = llt_denom.solve(Eigen::MatrixXd::Identity(p, p));
    double sigma = (Y - X * B).squaredNorm() / (Y.rows() - p);
    Eigen::VectorXd SE = sqrt(sigma) * (INV.selfadjointView<Eigen::Upper>() * Xt_PZ_X * INV.selfadjointView<Eigen::Upper>()).diagonal().cwiseSqrt();

    // Eigen::MatrixXd output(2, B.rows());
    // output << B.transpose(), SE.transpose();
    // return output;
    return List::create(Named("BETA") = B,
                        Named("SE") = SE,
                        Named("R2") = 1 - sigma * (Y.rows() - p) / (Y.rows() - 1));
  }
