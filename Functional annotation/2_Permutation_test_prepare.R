##########################3
# Prepare data for spatial permutation.

library(neurobase)

# read template nii
atlas <- readNIfTI("./neuromorphometrics.nii")

mni <- readNIfTI("/template_1.5mm/Template_1_IXI555_MNI152_GS.nii")
mni <- copyNIfTIHeader(mni, mni[,,,1])
mni@bitpix=32
mni@datatype=16
mni@data_type="16"

 # EUR only baseline VBM results
load(file="./EUR_voxel_GM_output.Rdata")

PVAL_base <- mni
PVAL_base[voxel_list] <- output$F
PVAL_base[-voxel_list] <- NA

writeNIfTI(PVAL_base, "./MAP_TO_PERM")
# resample
MAP <- fslr::flirt("./MAP_TO_PERM.nii.gz",  "./response_execution_map.nii", dof=12, omat="./omat", opts="-applyxfm -usesqform  -interp nearestneighbour")

########################33#########

coord <- which(!is.na(MAP), arr.ind = TRUE)
colnames(coord) <- c("x","y","z")

fwrite(coord, "./COORD_TO_PERM.txt", sep="\t", col.names=FALSE)
write.table(MAP[coord], "./SES_VAL_TO_PERM.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

#######################################
#######################################
#######################################





#######################################
# Repate for VBM of PGI and VBM of SES PCs with PGI controlled for. 

 # VBM of PGI
load(file="./PGS_voxel_GM_output.Rdata")

PVAL_base <- mni
PVAL_base[voxel_list] <- output$F
PVAL_base[-voxel_list] <- NA

writeNIfTI(PVAL_base, "./PGS_MAP_TO_PERM")
MAP <- fslr::flirt("./PGS_MAP_TO_PERM.nii.gz",  "./response_execution_map.nii", dof=12, omat="../../TEMP/omat", opts="-applyxfm -usesqform  -interp nearestneighbour")

coord <- which(!is.na(MAP), arr.ind = TRUE)
colnames(coord) <- c("x","y","z")
write.table(MAP[coord], "./PGS_VAL_TO_PERM.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


# VBM of SES PCs with PGI controlled for. 
load(file="../../OUTPUT/SESPC_PGS_voxel_GM_output.Rdata")

PVAL_base <- mni
PVAL_base[voxel_list] <- output$F
PVAL_base[-voxel_list] <- NA

writeNIfTI(PVAL_base, "./COND_PGS_MAP_TO_PERM")
MAP <- fslr::flirt("./COND_PGS_MAP_TO_PERM.nii.gz",  "./response_execution_map.nii", dof=12, omat="./omat", opts="-applyxfm -usesqform  -interp nearestneighbour")

coord <- which(!is.na(MAP), arr.ind = TRUE)
colnames(coord) <- c("x","y","z")
write.table(MAP[coord], "./COND_PGS_VAL_TO_PERM.txt", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)


