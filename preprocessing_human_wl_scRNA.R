#library(BiocManager)
#library(remotes)
#library(devtools)
#devtools::install_github('satijalab/seurat-wrappers')
#library(R.utils)
#library(SeuratWrappers)
install.packages("rlang", dependencies = TRUE)

library(monocle3)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(patchwork)
library(sctransform)
library(glmGamPoi)


p1_t0 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7142\\filtered_feature_bc_matrix")
p1_t0 <- CreateSeuratObject(counts = p1_t0, project = "p1")
p1_t1 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7143\\filtered_feature_bc_matrix")
p1_t1 <- CreateSeuratObject(counts = p1_t1, project = "p1")
p1_t2 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7144\\filtered_feature_bc_matrix")
p1_t2 <- CreateSeuratObject(counts = p1_t2, project = "p1")
p2_t0 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7145\\filtered_feature_bc_matrix")
p2_t0 <- CreateSeuratObject(counts = p2_t0, project = "p2")
p2_t1 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7146\\filtered_feature_bc_matrix")
p2_t1 <- CreateSeuratObject(counts = p2_t1, project = "p2")
p2_t2 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7147\\filtered_feature_bc_matrix")
p2_t2 <- CreateSeuratObject(counts = p2_t2, project = "p2")
p3_t0 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7148\\filtered_feature_bc_matrix")
p3_t0 <- CreateSeuratObject(counts = p3_t0, project = "p3")
p3_t1 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7149\\filtered_feature_bc_matrix")
p3_t1 <- CreateSeuratObject(counts = p3_t1, project = "p3")
p3_t2 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7150\\filtered_feature_bc_matrix")
p3_t2 <- CreateSeuratObject(counts = p3_t2, project = "p3")
p4_t0 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7151\\filtered_feature_bc_matrix")
p4_t0 <- CreateSeuratObject(counts = p4_t0, project = "p4")
p4_t1 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7152\\filtered_feature_bc_matrix")
p4_t1 <- CreateSeuratObject(counts = p4_t1, project = "p4")
p4_t2 <- Read10X(data.dir ="C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\D19-7153\\filtered_feature_bc_matrix")
p4_t2 <- CreateSeuratObject(counts = p4_t2, project = "p4")

adata <- merge(p1_t0, y = c(p1_t1,p1_t2,p2_t0,p2_t1,p2_t2,p3_t0,p3_t1,p3_t2,p4_t0,p4_t1,p4_t2), 
                       add.cell.ids = c('p1_t0','p1_t1','p1_t2','p2_t0','p2_t1','p2_t2','p3_t0','p3_t1','p3_t2','p4_t0','p4_t1','p4_t2'), 
                       project = "human_wl")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
adata[["percent.mt"]] <- PercentageFeatureSet(object = adata, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(adata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

adata <- SCTransform(adata, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE, conserve.memory=TRUE)

# These are now standard steps in the Seurat workflow for visualization and clustering
adata <- RunPCA(adata, verbose = FALSE)
adata <- RunUMAP(adata, dims = 1:30, verbose = FALSE)

adata <- FindNeighbors(adata, dims = 1:30, verbose = FALSE)
adata <- FindClusters(adata, verbose = FALSE)
DimPlot(adata, label = TRUE) + NoLegend()

saveRDS(adata, file = "C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Human Weight Loss Cohort\\WL_PATIENT_10X_Human\\human_wl_adata.rds")