# Testing differences in SES-GMV associations with and without PGI as a control variable (Supplemenatry Section 4.3.2 )

##########################################################################################
UKB <- readRDS("./UKB_Pheno.Rds")
UKB <- UKB[!is.na(SES_GSEM_PGS)]

##########################################################################################
# construct SES PC
# combine neighborhood srores
    score_list <- c("inc_score", "emp_score", "edu_score")
    std_list <- paste0(score_list, "_std")
    UKB[, (std_list) := lapply(.SD, function(x){(x-mean(x))/sd(x)}), .SDcols=score_list ]
    
    cov <- cov(UKB[, .SD, .SDcols=std_list])

    one=matrix(rep(1, length(std_list)), nrow=length(std_list))
    denom=as.numeric( (t(one) %*% solve(cov)  %*% one) )
    W=(t(one) %*% solve(cov)) / denom

    UKB[, score := apply(.SD, 1, function(x){W %*% x}), .SDcols=std_list ]
    UKB[, score := score / sd(score)]


library(PCAmixdata)
PC <- PCAmix(X.quanti=UKB[sample_SES==1, .(score, log_soc_mean, log_MSOA_y)], X.quali=UKB[sample_SES==1, .(edu, HHinc, house, SOC_3d)], rename.level=T, graph=F, ndim=25)

ndim=2
UKB[sample_SES==1, (paste0("SES_", 1:ndim)) := data.table(PC$scores[,1:ndim]) ] 

# change sign
UKB[, SES_1 := - SES_1]

# standardize SES measrue
UKB[, SES_1:= (SES_1 - mean(SES_1)) / sd(SES_1)] 
UKB[, SES_2:= (SES_2 - mean(SES_2)) / sd(SES_2)] 

#########################################################################################
# standardize #
PGS_list <- grep("PGS", names(UKB), value=T)
UKB[, (PGS_list) := lapply(.SD, function(x) (x-mean(x))/sd(x)) , .SDcols= PGS_list]

##########################################################################################
pc <- paste0("PC", 1:40)
pc <- paste(pc, collapse="+")
SES <- paste0("SES_", 1:2)
SES <- paste0(SES, collapse="+")
fm <- paste0("SES_GSEM_PGS ~", SES , "+ ( (age + I(age^2) + I(age^3))*male + TIV + date_sp_1+date_sp_2+date_sp_3 + time)*site + genop  +", pc)

lm <- s(lm(as.formula(fm), UKB))

b <- matrix(lm$coefficients[2:3,1])
b_cov <- vcov(lm)[2:3, 2:3]
btb  <- b %*% t(b)
##########################################################################################

# Load VBM output on SES
load(paste0("./EUR_voxel_GM_output.Rdata"))
T_beta_1 = output$beta_1 
T_beta_2 = output$beta_2
sig = output$sig_5
p_rsq_0 = output$p_rsq
SSR_diff = output$SSR_r - output$SSR

# Load VBM output on SES + PGI
load(paste0("./SESPC_PGS_voxel_GM_output.Rdata"))

get_stat <- function(a, se){
V = se^2*btb + a^2*b_cov
a^2 * t(b) %*% solve(V) %*% b
}

output[, sig_AFTER := sig_5]
output[, sig_BEFORE := sig]

output[, T_beta_1 := T_beta_1]
output[, T_beta_2 := T_beta_2]

output[, diff_1 := beta_3*b[1]]
output[, diff_2 := beta_3*b[2]]

output[, diff_1_pc := diff_1/T_beta_1]
output[, diff_2_pc := diff_2/T_beta_2]

output[, p_rsq_diff := (SSR_diff - (SSR_r - SSR) )/ SSR_diff ] # % due to including PGI

output[, eq_chi:= mapply(get_stat, output$beta_3, output$se_3)]
output[, eq_pval := pchisq(eq_chi, 2, lower.tail=FALSE)]






