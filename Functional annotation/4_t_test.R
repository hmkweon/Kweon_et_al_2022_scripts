library(neurobase)
load("./ADJ_PVAL_CUT.Rdata") # FWE p-value cutoff

NCORE = 100

############################################33

LIST <- fread("./cognitive_atlas_concept_category.csv")
# exclude problematic ones
to_drop <- as.character(t(fread("./drop_list.csv", header=FALSE)))
to_drop <- to_drop[- which(to_drop %in% c("of","to","down","under"))]
to_drop <- paste0(to_drop, collapse="|")
to_drop <- LIST[, grep(to_drop, concept, value=TRUE) ]
# 12

LIST <- LIST[!concept %in% to_drop]

STAT = lapply(1:nrow(LIST), function(x){
name <- paste0(LIST[x, stringr::str_split(concept, " ")][[1]], collapse="_")
temp <- RNifti::readNifti(paste0("./", name, "_map.nii"))
temp <- temp^2
data.table(MEAN = mean(temp, na.rm=TRUE), MEDIAN = median(temp, na.rm=TRUE), MAX=max(temp, na.rm=TRUE), MIN=min(temp, na.rm=TRUE), PERC_SIG = sum(pchisq(temp, df=1, lower.tail=FALSE) < 0.05) / length(temp) )
})

LIST = cbind(LIST, rbindlist(STAT))
LIST = LIST[MEAN != 0]
############################################33

# collect maps
get_maps <- function(x){
name <- paste0(LIST[x, stringr::str_split(concept, " ")][[1]], collapse="_")
temp <- RNifti::readNifti(paste0("./NEUROQUERY/", name, "_map.nii"))

return(temp[COORD]^2)
}

COORD = as.matrix(fread("./COORD_TO_PERM.txt"))

flexiblas_set_num_threads(1)
MAPS <- pbmcapply::pbmcmapply( get_maps, 1:nrow(LIST), mc.cores=NCORE, ignore.interactive=T)
flexiblas_set_num_threads(NCORE)

MAPS <- data.table(MAPS)

RM_IND = NULL
for (i in 1:nrow(MAPS)){
if (all(as.numeric(MAPS[i, ])==0)) RM_IND = c(RM_IND, i)
}

MAPS = MAPS[!RM_IND]

################################################
t_func <- function(x, NULL_MAP){
temp = MAPS[[x]]
EST = coef(s(lm(temp ~ SIG[ INCL_IND])))[2,3]

PERM = sapply(1:10000, function(i){
coef(s(lm(temp ~ ifelse(pf(NULL_MAP[i, INCL_IND], 2, 20717, lower.tail=FALSE) < 0.01, 1, 0) )))[2,3] } )

PVAL = sum(PERM >=  EST) / length(PERM)
return(list(EST = EST, PERM = PERM, PVAL=PVAL))
gc(verbose=FALSE)
}
##################################3
# Bseline VBM (SES PC)
NULL_MAP <- as.matrix(fread("./NULL_SES_MAPS.txt"))

SIG = fread("./SES_VAL_TO_PERM.txt")[[1]]
SIG = ifelse(pf(SIG, 2, 20717, lower.tail=FALSE) < 0.01, 1, 0)

INCL_IND = 1:length(SIG)
INCL_IND = INCL_IND[!INCL_IND %in% RM_IND]

gc(verbose=FALSE)

flexiblas_set_num_threads(1)
RES <- pbmcapply::pbmclapply(1:nrow(LIST), t_func, NULL_MAP=NULL_MAP, mc.cores=NCORE, ignore.interactive=T)
flexiblas_set_num_threads(NCORE)

PERM_MAX = sapply(1:10000, function(i) max(sapply(1:nrow(LIST), function(x) RES[[x]]$PERM[i])))

LIST$SES_T <- sapply(1:nrow(LIST), function(i) RES[[i]]$EST)
LIST$SES_T_pval <- sapply(1:nrow(LIST), function(i) RES[[i]]$PVAL)
LIST[order(-abs(SES_T))][1:50]
LIST[order(SES_T_pval)][1:50]

LIST[SES_T_pval<0.05]
LIST[p.adjust(SES_T_pval, "fdr")<0.05]

LIST$SES_T_pval_perm  <- sapply(1:nrow(LIST), function(j) sum(sapply(1:10000, function(i) LIST$SES_T[j] <= PERM_MAX[i]) ) / length(PERM_MAX)  )

##############################
# VBM of PGI
NULL_MAP <- as.matrix(fread("./NULL_PGS_MAPS.txt"))

SIG = fread("./PGS_VAL_TO_PERM.txt")[[1]]
SIG = ifelse(pf(SIG, 2, 20717, lower.tail=FALSE) < 0.001, 1, 0)

gc(verbose=FALSE)

flexiblas_set_num_threads(1)
RES <- pbmcapply::pbmclapply(1:nrow(LIST), t_func, NULL_MAP=NULL_MAP, mc.cores=NCORE, ignore.interactive=T)
flexiblas_set_num_threads(NCORE)

LIST$PGS_T <- sapply(1:nrow(LIST), function(i) RES[[i]]$EST)
LIST$PGS_T_pval <- sapply(1:nrow(LIST), function(i) RES[[i]]$PVAL)
LIST[order(-abs(PGS_T))][1:50]
LIST[order(PGS_T_pval)][1:50]

LIST[PGS_T_pval<0.05]
LIST[p.adjust(PGS_T_pval, "fdr")<0.05]

PERM_MAX = sapply(1:10000, function(i) max(sapply(1:nrow(LIST), function(x) RES[[x]]$PERM[i])))
LIST$PGS_T_pval_perm <- sapply(1:nrow(LIST), function(j) sum(sapply(1:10000, function(i) LIST$PGS_T[j] <= PERM_MAX[i]) ) / length(PERM_MAX)  )

gc(verbose=FALSE)

##############################
# VBM of SES PC with PGI control
NULL_MAP <- as.matrix(fread("./NULL_COND_PGS_MAPS.txt"))
SIG = fread("./COND_PGS_VAL_TO_PERM.txt")[[1]]
SIG = ifelse(SIG > 13.03721, 1, 0)

gc(verbose=FALSE)

flexiblas_set_num_threads(1)
RES <- pbmcapply::pbmclapply(1:nrow(LIST), t_func, NULL_MAP=NULL_MAP, mc.cores=NCORE, ignore.interactive=T)
flexiblas_set_num_threads(NCORE)

LIST$COND_PGS_T <- sapply(1:nrow(LIST), function(i) RES[[i]]$EST)
LIST$COND_PGS_T_pval <- sapply(1:nrow(LIST), function(i) RES[[i]]$PVAL)
LIST[order(-abs(COND_PGS_T))][1:50]
LIST[order(COND_PGS_T_pval)][1:50]

LIST[COND_PGS_T_pval<0.05]
LIST[p.adjust(COND_PGS_T_pval, "fdr")<0.05]

PERM_MAX = sapply(1:10000, function(i) max(sapply(1:nrow(LIST), function(x) RES[[x]]$PERM[i])))
LIST$COND_PGS_T_pval_perm <- sapply(1:nrow(LIST), function(j) sum(sapply(1:10000, function(i) LIST$COND_PGS_T[j] <= PERM_MAX[i]) ) / length(PERM_MAX)  )

gc(verbose=FALSE)

fwrite(LIST, "./NEUROQUERY_T_output.csv", sep="\t")
