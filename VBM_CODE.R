#### R Script to run the baseline VBM 
#### Written by Hyeokmoon Kweon (h.kweon@vu.nl)

## read the phenotype data
UKB <- readRDS("./UKB_Pheno.Rds")

##############################################################################################
##############################################################################################
# Baseline VBM
##############################################################################################
##############################################################################################

## prepare the models
pc <- paste0("PC", 1:40)
pc <- paste(pc, collapse="+")
SES <- paste0("SES_", 1:2)
SES <- paste0(SES, collapse="+")
fm <- paste0("~", SES , " + ( (age + I(age^2) + I(age^3))*male + TIV  + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)
fm_r <- paste0("~  ( (age + I(age^2) + I(age^3))*male + TIV  + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)

X <- model.matrix( as.formula(fm), UKB )
X_r <- model.matrix( as.formula(fm_r), UKB )

##QR decomposition
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)
XX <- chol2inv(R)

varX <- diag(XX)[grepl("SES", colnames(X))]
DF <- nrow(UKB) - QR$rank
varX_div_DF <- varX / DF 

QR_r <- qr(X_r)
Q_r <- qr.Q(QR_r)


######## VBM run in batch (START) #######

# batch number (1 as example)
i = 1

library(CppFunc) # personal C++ estimation codes https://github.com/hmkweon/CppFunc
library(pbmcapply)


# Read grey matter volume data
data <- read_fst("")

# get chunks to run parallel over
chunk_s <- seq(1, ncol(data), by=500)
chunk_e <- c(chunk_s[2:length(chunk_s)]-1, ncol(data))

# function to run in parallel
run_reg <- function(x){
reg_F4(R, Q, Q_r, data[,chunk_s[x]:chunk_e[x]] , sample_ind, length(varX_div_DF))
}

# Run
result <- pbmclapply(1:length(chunk_e), run_reg,  mc.cores=20,  ignore.interactive=T)

# collect results
output <- matrix(NA, length(varX_div_DF)+2, ncol(data))
names(output) <- names(data)

for (j in 1:length(chunk_e)){
output[, chunk_s[j]:chunk_e[j]] <- result[[j]]
}

output <- data.table(output)
file_out <- paste0("../OUTPUT/GM_", i, ".fst")
write_fst(output, file_out, compress=0)

######## VBM run in batch (END) #######

## Collect the results from each batch and compute test statistics and effect sizes
ncol=504426
s_col <- seq(1, ncol, by=30000)
e_col <- c(s_col[2:length(s_col)]-1, ncol)

output <- matrix(NA, nrow=ncol, ncol=length(varX_div_DF)+2) 
for (i in 1:length(s_col)){
file <- paste0("../../OUTPUT/GM_", i, ".fst")
output[s_col[i]:e_col[i],] <- t(as.matrix(read_fst(file)))
}

output <- data.table(output)
setnames(output, names(output), c(paste0("beta_", 1:length(varX_div_DF)),  "SSR", "SSR_r") )

for (i in 1:length(varX_div_DF)){
b <- paste0("beta_", i)
se <- paste0("se_", i)
t <- paste0("t_", i)
output[, (se) := sqrt(SSR * varX_div_DF[i])  ]
output[, (t) := get(b)/get(se) ]
output[, paste0("pval_", i) := 2*pt(-abs(get(t)), df=DF)]
output[, paste0("p_rsq_",i) := get(t)^2 / ( get(t)^2 + DF) ] # partial R2
}
output[, p_rsq := (SSR_r-SSR)/SSR_r] 
output[, F := ((SSR_r-SSR)/2) / (SSR/DF)  ]
output[, F_pval := pf(F, 2, DF, lower.tail=FALSE)]

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################




##############################################################################################
##############################################################################################
# VBM on SES PC controlling for PGI (=PGS) with GIV
##############################################################################################
##############################################################################################

# Sample with European ancestry only
UKB <- UKB[!is.na(SES_GSEM_PGS)]


############################################################################33
## prepare the models
fm_x <- paste0("~ SES_1 + SES_2 + SES_GSEM_PGS_sub1 + ( (age + I(age^2) + I(age^3))*male + TIV + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)
fm_z <- paste0("~ SES_1 + SES_2 + SES_GSEM_PGS_sub2 + ( (age + I(age^2) + I(age^3))*male + TIV + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)
fm_r <- paste0("~ + SES_GSEM_PGS_sub1 + ( (age + I(age^2) + I(age^3))*male + TIV + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)

X <- model.matrix( as.formula(fm_x), UKB )
Z <- model.matrix( as.formula(fm_z), UKB )

##QR decomposition
Q_Z <- qr.Q(qr(Z))
Xh <- (Q_Z %*% crossprod(Q_Z, X))

QR <- qr(Xh)
Q <- qr.Q(QR)
R <- qr.R(QR)
XX <- chol2inv(R)

varX <- diag(XX)[grepl("PGS|SES", colnames(X))]
DF <- nrow(UKB) - QR$rank
varX_div_DF <- varX / DF 

X_r <- model.matrix( as.formula(fm_r), UKB )
Xh_r <- (Q_Z %*% crossprod(Q_Z, X_r))
Q_r <- qr.Q(qr(Xh_r))

######## VBM run in batch (START) #######

# batch number (1 as example)
i = 1

library(CppFunc) # personal C++ estimation codes https://github.com/hmkweon/CppFunc
library(pbmcapply)

# Read grey matter volume data
data <- read_fst("")

# get chunks to run parallel over
chunk_s <- seq(1, ncol(data), by=500)
chunk_e <- c(chunk_s[2:length(chunk_s)]-1, ncol(data))

# function to run in parallel
run_reg <- function(x){
IV_F(X, R, Q, Q_r, data[[x]] , sample_ind, length(varX_div_DF))
}

# Run
result <- pbmclapply(1:length(chunk_e), run_reg,  mc.cores=20,  ignore.interactive=T)

# collect results
output <- matrix(NA, length(varX_div_DF)+2, ncol(data))
names(output) <- names(data)

for (j in 1:length(chunk_e)){
output[, chunk_s[j]:chunk_e[j]] <- result[[j]]
}

output <- data.table(output)
file_out <- paste0("../OUTPUT/GM_", i, "_GIV.fst")
write_fst(output, file_out, compress=0)

######## VBM run in batch (END) #######

## Collect the results from each batch and compute test statistics and effect sizes
ncol=504426
s_col <- seq(1, ncol, by=30000)
e_col <- c(s_col[2:length(s_col)]-1, ncol)

output <- matrix(NA, nrow=ncol, ncol=length(varX_div_DF)+2) 
for (i in 1:length(s_col)){
file <- paste0("../../OUTPUT/GM_", i, ".fst")
output[s_col[i]:e_col[i],] <- t(as.matrix(read_fst(file)))
}


output <- data.table(output)
setnames(output, names(output), c(paste0("beta_", 1:length(varX_div_DF)), "SSR_2sls", "SSR", "SSR_r") )

for (i in 1:length(varX_div_DF)){
b <- paste0("beta_", i)
se <- paste0("se_", i)
t <- paste0("t_", i)
output[, (se) := sqrt(SSR_2sls * varX_div_DF[i])  ]
output[, (t) := get(b)/get(se) ]
output[, paste0("pval_", i) := 2*pt(-abs(get(t)), df=DF)]
output[, paste0("d_rsq_",i) := (get(t)^2 * SSR_2sls) / ((nrow(Q)-1)*DF) ]   # fine for PGS, too. 
output[, paste0("p_rsq_",i) := get(t)^2 / ( get(t)^2 + DF) ]   
}
#
output[, p_rsq := (SSR_r-SSR)/SSR_r] 
output[, F := ((SSR_r-SSR)/ 2) / (SSR_2sls/DF)   ]
output[, F_pval := pf(F,  2, DF, lower.tail=FALSE)]

##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################
