## Functions for cluster-level analysis
run_reg <- function(fm, fm_r, var, Y, data=UKB){
library(CppFunc)
X <- model.matrix( as.formula(fm), data )
X_r <- model.matrix( as.formula(fm_r), data )

X <- cbind(X[, c(1,grep(var, colnames(X)))], X[, -c(1,grep(var, colnames(X)))])

QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)
XX <- chol2inv(R)

varX <- diag(XX)[grepl(var, colnames(X))]
DF <- nrow(data) - QR$rank
varX_div_DF <- varX / DF 

QR_r <- qr(X_r)
Q_r <- qr.Q(QR_r)

res <- reg_F_ROI(R, Q, Q_r, Y, length(varX_div_DF))

out <- data.table(t(res))
setnames(out, names(out), c(paste0("beta_", 1:length(varX_div_DF)), "SSR", "SSR_r"))

for (i in 1:length(varX_div_DF)){
b <- paste0("beta_", i)
se <- paste0("se_", i)
t <- paste0("t_", i)
out[, (se) := sqrt(SSR * varX_div_DF[i])  ]
out[, (t) := get(b)/get(se) ]
out[, paste0("pval_", i) := 2*pt(-abs(get(t)), df=DF)]
out[, paste0("d_rsq_",i) := get(t)^2 * (SSR/(nrow(Q)-1)) / DF]
}
#
out[, d_rsq := (SSR_r-SSR)/(nrow(Q)-1)] # length(Y)-1
out[, p_rsq := (SSR_r-SSR)/SSR_r] # length(Y)-1

if(length(varX_div_DF)>1){
out[, F := ((SSR_r-SSR)/2) / (SSR/DF)   ]
out[, F_pval := pf(F, 2, DF, lower.tail=FALSE)]
} else if(length(varX_div_DF)==1){
out[, F := ((SSR_r-SSR)/1) / (SSR/DF)   ]
out[, F_pval := pf(F, 1, DF, lower.tail=FALSE)]
}

return(out)
}



run_iv <- function(fm_x, fm_z, fm_r, Y){

X <- model.matrix( as.formula(fm_x), UKB )
Z <- model.matrix( as.formula(fm_z), UKB )
Q_Z <- qr.Q(qr(Z))

Xh <- (Q_Z %*% crossprod(Q_Z, X))

QR <- qr(Xh)
Q <- qr.Q(QR)
R <- qr.R(QR)
XX <- chol2inv(R)

varX <- diag(XX)[grepl("SES", colnames(X))]
DF <- nrow(UKB) - QR$rank
varX_div_DF <- varX / DF 

X_r <- model.matrix( as.formula(fm_r), UKB )
Xh_r <- (Q_Z %*% crossprod(Q_Z, X_r))
Q_r <- qr.Q(qr(Xh_r))

Y <- scale(Y)

B <- qr.coef(QR, Y)
SSR_2sls <- colSums((Y - X%*%B)^2)
SSR <- colSums((qr.resid(QR, Y))^2)
SSR_r <- colSums((qr.resid(qr(Xh_r), Y))^2)

out <- data.table(t(rbind(B[2:4,], SSR_2sls, SSR, SSR_r)  ))
setnames(out, names(out), c(paste0("beta_", 1:length(varX_div_DF)), "SSR_2sls", "SSR", "SSR_r") )

for (i in 1:length(varX_div_DF)){
b <- paste0("beta_", i)
se <- paste0("se_", i)
t <- paste0("t_", i)

out[, (se) := sqrt(SSR_2sls * varX_div_DF[i])  ]
out[, (t) := get(b)/get(se) ]
out[, paste0("pval_", i) := 2*pt(-abs(get(t)), df=DF)]
out[, paste0("d_rsq_",i) := (get(t)^2 * SSR_2sls) / ((nrow(Q)-1)*DF) ]   # fine for PGS, too. 
out[, paste0("p_rsq_",i) := get(t)^2 / ( get(t)^2 + DF) ]   

}
#
out[, d_rsq := (SSR_r-SSR)/(nrow(Q)-1)] # length(Y)-1
out[, p_rsq := (SSR_r-SSR)/SSR_r] # length(Y)-1

out[, F := ((SSR_r-SSR)/ (2)) / (SSR_2sls/DF)   ]
out[, F_pval := pf(F,  (2), DF, lower.tail=FALSE)]

return(out)
}
