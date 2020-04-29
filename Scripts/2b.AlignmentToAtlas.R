# This part still under Seurat v2, like the initial analysis.
detach(package:Seurat, unload=TRUE) # These three lines are just here so that this part of the analysis can be re-run repeatedly, without re-running the part 1.Load_QC_CCA. 
library(Seurat,lib.loc = path_Seuratv2)
rm(SOH,SOH2)

# Create a Seurat Object for Harmony merging
SOH <- MergeSeurat(object1 = SO, object2 = SO.atlas)
unique(SOH@meta.data$group)
unique(SOH@meta.data$orig.ident)
unique(SOH@ident)
unique(SOH@meta.data$cell.type)

# Take the intersection of variable genes of the various datasets as the variable genes for the alignment:
Reduce(intersect,list(SO.C1@var.genes,SO.C2@var.genes,SO.C3@var.genes,SO.C4@var.genes,SO.atlas@var.genes)) -> var.genes.intersect
SOH@var.genes <- var.genes.intersect
length(SOH@var.genes)

# Compute PCA of the merged datasets:
SOH <- ScaleData(object = SOH)
SOH <- RunPCA(SOH, pcs.compute = 40)
PCElbowPlot(object = SOH, num.pc = 40)
dims.use <- 1:30

# From here on, switch to Seurat v3 (for compatibility with Harmony).
detach(package:Seurat, unload=TRUE)
library(Seurat)
library(harmony)
SOH2 <- UpdateSeuratObject(SOH)

tmp <- as.character(SOH2@meta.data$orig.ident)
tmp[tmp %in% c("B2","B3","B4","B5","B6","B7","B8","B9","B10")] <- "in vivo"
tmp[tmp %in% c("C1","C2","C3","C4")] <- "in vitro"
tmp <- factor(tmp)
unique(tmp)
SOH2@meta.data$dataset <- tmp
Idents(SOH2) <- tmp

SOH2 <- RunHarmony(object = SOH2, group.by.vars = c("dataset"), dims.use = dims.use, plot_convergence = T)
SOH2 <- RunUMAP(object = SOH2, reduction = "harmony", dims = dims.use)
DimPlot(object = SOH2, reduction = "umap")
SOH2 <- RunHarmony(object = SOH2, group.by.vars = c("orig.ident"), dims.use = dims.use, plot_convergence = T,reduction="harmony",reduction.save="harmony2")
#SOH2 <- RunUMAP(object = SOH2, reduction = "harmony2", dims = dims.use, reduction.name = "umap2",n.neighbors = 150, min.dist = 1, spread = 1,seed.use = 173)
SOH2 <- RunUMAP(object = SOH2, reduction = "harmony2", dims = dims.use, reduction.name = "umap2",n.neighbors = 150, min.dist = 1, spread = 1,seed.use = 700)
DimPlot(object = SOH2, reduction = "umap2",group.by = "cell.type",pt.size=1) # Export 12x8
DimPlot(object = SOH2, reduction = "umap2",group.by = "dataset",pt.size=1)
DimPlot(object = SOH2, reduction = "umap2",group.by = "orig.ident",pt.size=1)
#SOH2@reductions$umap <- NULL # Might need to run this line of code if running into a bug of Seurat using umap instead of umap2 in the feature plot despite of the explicit request for umap2.
FeaturePlot(object = SOH2, features = c("Lgr5","Lyz1","Tff3","Fabp1","Chga","Cdk1","Marcksl1","Ctgf","Anxa5"),reduction = "umap2", ncol = 3, dims = c(1,2),sort.cell = T)
FeaturePlot(object = SOH2,features = tmp_markers,reduction = "umap2", ncol = 3, dims = c(1,2),sort.cell = T)

# If one wants to go back to Seurat v2 to continue the main analysis:
detach(package:Seurat, unload=TRUE)
library(Seurat,lib.loc = path_Seuratv2)
