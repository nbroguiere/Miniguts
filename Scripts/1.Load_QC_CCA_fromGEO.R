##### Set paths #####
setwd("C:/Users/nicki/Documents/PostDoc_LSCB/19-03-20 Datasets Mike/CleanedUpCodeForUpload/Scripts")
path_Seuratv2 <- '../../Seuratv2' # Analysis was done in Seurat v2, and Seurat v3 is not backwards compatible so need a separate install of the old Seurat to keep running this code
path.data <- '../Resources/raw_counts.tsv' # Download and unzip GSE148366_raw_counts.tsv.gz
path.data.atlas <- 'C:/Users/nicki/Documents/PostDoc_LSCB/19-04-11 Single cell survey of mSI epithelium/GSE92332_atlas_UMIcounts.txt' # Haber 2017, in vivo atlas, download from GEO 

##### Library import #####
library(rPython)
library(Matrix)
library(Seurat,lib.loc = path_Seuratv2)
library(dplyr)
library(scater)        # SingleCellExperiment class
library(scRNAseq)
library(DropletUtils)  # read10xCounts function
library(matrixStats)
library(magrittr)
library(ggplot2)
#library(Rmagic)
#library(biomaRt)
#library(SCENIC)
#library(AUCell)
#library(Rtsne)
library(reticulate)
library(limma)
library(ggeffects)
register(SerialParam())
theme_set(theme_light())

##### Load raw data #####
raw.data <- as.matrix(read.csv(file = path.data, header = T, sep = "\t"))
atlas <- as.matrix(read.table(file = path.data.atlas, sep='\t'))

##### QC, normalize #####
SO.tmp <- CreateSeuratObject(raw.data = raw.data, min.cells = 4, min.genes = 300, project = "OrganoidsAndMiniguts")
mito.genes       <- grep(pattern = "^mt-", x = rownames(x = SO.tmp@data), value = TRUE)
percent.mito     <- Matrix::colSums(SO.tmp@raw.data[mito.genes, ])/Matrix::colSums(SO.tmp@raw.data)
percent.mito.log <- log10(percent.mito*100)
nGene.log        <- log10(SO.tmp@meta.data$nGene)
names(nGene.log) <- rownames(SO.tmp@meta.data)
nUMI.log         <- log10(SO.tmp@meta.data$nUMI)
names(nUMI.log) <- rownames(SO.tmp@meta.data)
hist(SO.tmp@meta.data$nGene,100,col="blue")
hist(percent.mito,100,col="blue")
hist(SO.tmp@meta.data$nUMI,200,col="blue")  # 
SO.tmp <- AddMetaData(object = SO.tmp, metadata = percent.mito, col.name = "percent.mito")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = percent.mito.log, col.name = "percent.mito.log")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = nGene.log, col.name = "nGene.log")
SO.tmp <- AddMetaData(object = SO.tmp, metadata = nUMI.log, col.name = "nUMI.log")
ggplot(SO.tmp@meta.data) + geom_point(aes(x=nGene, y=percent.mito, color=nUMI.log))+ scale_x_continuous(trans='log10')+scale_y_continuous(trans='log10')
ggplot(SO.tmp@meta.data) + geom_point(aes(x=nGene, y=percent.mito, color=nUMI.log))
ggplot(SO.tmp@meta.data) + geom_point(aes(x=nGene.log, y=nUMI.log, color=nUMI.log))
SO.tmp <- FilterCells(object = SO.tmp, subset.names = c("nGene", "percent.mito"), low.thresholds = c(2500, 0.04), high.thresholds = c(5000, 0.15))
SO.tmp <- NormalizeData(object = SO.tmp, normalization.method = "LogNormalize", scale.factor = 10000)

SO.atlas <- CreateSeuratObject(raw.data=atlas,min.cells=1,min.genes=300,project="atlas")
nGene.log        <- log10(SO.atlas@meta.data$nGene)
names(nGene.log) <- rownames(SO.atlas@meta.data)
nUMI.log         <- log10(SO.atlas@meta.data$nUMI)
names(nUMI.log) <- rownames(SO.atlas@meta.data)
nUMI.nGene.logRatio <- nUMI.log-nGene.log
names(nUMI.nGene.logRatio) <- rownames(SO.atlas@meta.data)
hist(SO.atlas@meta.data$nGene,100,col="blue")  # Single live cells: 3200 - 4600
hist(SO.atlas@meta.data$nUMI,200,col="blue")  # 
SO.atlas <- AddMetaData(object = SO.atlas, metadata = nUMI.nGene.logRatio, col.name = "nUMI.nGene.logRatio")
SO.atlas <- AddMetaData(object = SO.atlas, metadata = nGene.log, col.name = "nGene.log")
SO.atlas <- AddMetaData(object = SO.atlas, metadata = nUMI.log, col.name = "nUMI.log")
hist(SO.atlas@meta.data$nUMI.nGene.logRatio,200,col="blue")  # 
ggplot(SO.atlas@meta.data) + geom_point(aes(x=nGene.log, y=nUMI.log, color=nUMI.log))
SO.atlas <- NormalizeData(object = SO.atlas, normalization.method = "LogNormalize", scale.factor = 10000)
SO.atlas <- SubsetData(SO.atlas, ident.use = paste0("B",2:10), subset.raw = F) # B1 to B10 are the various mice/datasets from Haber et al 2017. B1 has very strong technical differences, doesn't co-UMAP well with the rest without alignment, exclude for simplicity. 

##### Pick up the metadata from column names, especially cell type annotations #####
Dataset_CellType_UniqueTag <- strsplit2(SO.tmp@cell.names,split="_")
cell.type <- Dataset_CellType_UniqueTag[,2]
cell.type[cell.type=="G2.M"] <- "G2/M"
cell.type <- factor(cell.type, ordered = T, levels = c("ISC","G2/M","Enterocytes","TopEnterocytes","Microfold","Paneth","Goblet","Tuft","Enteroendocrine")) # ISC are the stem and progenitor, G2/M dividing cells, TopEnterocytes show analogies with Villus top enterocytes, Microfold are the RSC/Fetal like regenerative/M-like cells
names(cell.type) <- SO.tmp@cell.names
SO.tmp@meta.data$cell.type <- cell.type

##### Rename the datasets C1 for Organoid, C2 for MinigutYoung, C3 for MinigutOld, C4 for MinigutCrypto #####
tmp <- Dataset_CellType_UniqueTag[,1]
SO.tmp@meta.data$sample.name <- tmp
tmp[tmp=="Organoid"] <- "C1"
tmp[tmp=="MinigutYoung"] <- "C2"
tmp[tmp=="MinigutOld"] <- "C3"
tmp[tmp=="MinigutCrypto"] <- "C4"
tmp <- factor(tmp, levels=paste0("C",1:4))
names(tmp) <- SO.tmp@cell.names
SO.tmp@meta.data$orig.ident <- tmp
SO.tmp@ident <- tmp

##### Split up the 4 datasets #####
SO.C1 <- SubsetData(object = SO.tmp, ident.use = "C1", subset.raw = T)
SO.C2 <- SubsetData(object = SO.tmp, ident.use = "C2", subset.raw = T)
SO.C3 <- SubsetData(object = SO.tmp, ident.use = "C3", subset.raw = T)
SO.C4 <- SubsetData(object = SO.tmp, ident.use = "C4", subset.raw = T)

##### Merging C1 to C4 with CCA #####
genes.C1 <- rownames(SO.C1@data)
SO.C1 <- FindVariableGenes(object = SO.C1, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = -0.1)
genes.C2 <- rownames(SO.C2@data)
SO.C2 <- FindVariableGenes(object = SO.C2, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = -0.1)
genes.C3 <- rownames(SO.C3@data)
SO.C3 <- FindVariableGenes(object = SO.C3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = -0.1)
genes.C4 <- rownames(SO.C4@data)
SO.C4 <- FindVariableGenes(object = SO.C4, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = -0.1)

Reduce(intersect,list(SO.C1@var.genes,SO.C2@var.genes,SO.C3@var.genes,SO.C4@var.genes)) -> var.genes.intersect
Reduce(union,list(SO.C1@var.genes,SO.C2@var.genes,SO.C3@var.genes,SO.C4@var.genes)) -> var.genes.union

length(genes.C1)
length(SO.C1@var.genes)
length(genes.C2)
length(SO.C2@var.genes)
length(genes.C3)
length(SO.C3@var.genes)
length(genes.C4)
length(SO.C4@var.genes)
var.genes.intersect %>% length
var.genes.union %>% length

SO.C1@meta.data$group <- "C1"
SO.C2@meta.data$group <- "C2"
SO.C3@meta.data$group <- "C3"
SO.C4@meta.data$group <- "C4"

SO.C1 <- ScaleData(SO.C1,do.scale=TRUE,do.center=TRUE,check.for.norm=TRUE)
SO.C2 <- ScaleData(SO.C2,do.scale=TRUE,do.center=TRUE,check.for.norm=TRUE)
SO.C3 <- ScaleData(SO.C3,do.scale=TRUE,do.center=TRUE,check.for.norm=TRUE)
SO.C4 <- ScaleData(SO.C4,do.scale=TRUE,do.center=TRUE,check.for.norm=TRUE)

tmp <- list(SO.C1,SO.C2,SO.C3,SO.C4)
SO <- RunMultiCCA(object.list = tmp, genes.use =var.genes.intersect, add.cell.ids = c("C1","C2","C3","C4"), num.ccs = 30)
#MetageneBicorPlot(SO, grouping.var = "group", dims.eval = 1:30, display.progress = FALSE)
dims.align <- 1:30
SO <- AlignSubspace(SO, reduction.type = "cca", grouping.var = "group", dims.align = dims.align)

cells <- colnames(SO@data)
cells.C1 <- cells[grep("^C1_",cells)]
cells.C2 <- cells[grep("^C2_",cells)]
cells.C3 <- cells[grep("^C3_",cells)]
cells.C4 <- cells[grep("^C4_",cells)]
ident_backup_groups <- SO@ident
genes <- rownames(SO@data)

##### UMAP non-linear embedding #####
dims.use <- 1:11
SO <- RunUMAP(SO, dims.use=dims.use,reduction.use="cca.aligned",n_neighbors=10,min_dist=0.7,seed.use = 62)

# Plot datasets
SO@ident <- factor(SO@meta.data$sample.name)
names(SO@ident) <- SO@cell.names
DimPlot(SO,reduction.use = 'umap',pt.size=1)

# Plot cell types
SO@ident <- factor(SO@meta.data$cell.type)
names(SO@ident) <- SO@cell.names
DimPlot(SO,reduction.use = 'umap',pt.size=1)

##### Pick up a set of cells and markers or QC metadata for highlight in UMAP #####
tmp_markers <- c("Olfm4","Lgr5","Tff3","Fabp1","Mki67","Lyz1","Defa17","Chga","Krt19","Krt20","Dclk1","Anxa5")
tmp_markers <- c("nGene.log","nGene","percent.mito.log","percent.mito")
tmp_cells <- c(cells.C1, cells.C2, cells.C3, cells.C4)
FeaturePlot(object = SO, cells.use=tmp_cells, features.plot = intersect(tmp_markers,genes), cols.use = c("#EEEEEE", "#880EC9"), reduction.use = "umap",pt.size=1,no.axes = TRUE,no.legend = FALSE)

##### Gene set scoring, on merged in vitro datasets: Distal vs proximal enterocyte annotation and others e.g. microfold cells and defas #####
# Immature enterocyte markers from Haber 2017. Top enterocytes from Moor 2018.
top_genes <- readLines(con="../Resources/Moor_Villi_top_landmark_genes.txt")
bot_genes <- readLines(con="../Resources/Moor_Villi_bot_landmark_genes.txt")
M_genes <- readLines(con="./Resources/Haber_Mcell_landark_genes.txt")
top_genes <- intersect(top_genes,genes)
bot_genes <- intersect(bot_genes,genes)
M_genes <- intersect(M_genes,genes)
immature_genes <- c("Fgfbp1","Dmbt1","Pdss1","Prss32")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[top_genes,])), col.name = "top.genes")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[bot_genes,])), col.name = "bot.genes")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[immature_genes,])), col.name = "immature.genes")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[M_genes,])), col.name = "M.genes")
FeaturePlot(object = SO, features.plot = c("top.genes","bot.genes","immature.genes"), cols.use = c("#EEEEEE", "#880EC9"), reduction.use = "umap",pt.size=1)
top_genes_2 <- c("Ada","Ifrd1","Krt20","Pmp22","Serpinb1a")
bot_genes_2 <- c("Ccl25","Ckmt1","Cox4i1","Cox5a","Cox5b","Cox6b1","Cox6c","Cox7a1","Cox7a2","Cox7b","Gsta1","Fabp1","Gpx2","Gsta1","Lgals4","Lypd8","Minos1","Plac8","Reg3g","Sis","Spink4","Tma7","Txn1")
M_genes_2 <- c("Zmat3","Mmp15","Myadm","Anxa5","Marcksl1")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[top_genes_2,])), col.name = "top.genes.2")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[bot_genes_2,])), col.name = "bot.genes.2")
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[M_genes_2,])), col.name = "M.genes.2")
FeaturePlot(object = SO, features.plot = c(M_genes_2), cols.use = c("#EEEEEE", "#880EC9"), reduction.use = "umap",pt.size=1, no.legend = FALSE)

##### Cell cycle annotations from cell cycle genes published by the Regev lab, converted to mm: ##### 
cc.genes <- readLines(con="../Resources/Regev_cellcyclegenes_mm.txt")
s.genes <- cc.genes[1:42]
g2m.genes <- cc.genes[43:93]
SO@var.genes <- union(SO@var.genes,intersect(cc.genes,rownames(SO@data)))  ## Other useful function: setdiff(a,b) returns elements of a not in b. 
ident.backup <- SO@ident
SO <- CellCycleScoring(object = SO, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
SO <- AddMetaData(object = SO, metadata = ordered(SO@ident, levels = c("G1","S","G2M")), col.name = "cc.ident")
SO@ident <- ident.backup
# Plots on cell cycle: 
SO <- RunPCA(object = SO, pc.genes = c(s.genes, g2m.genes), do.print = FALSE,reduction.name="cc.pca")
DimPlot(object = SO,reduction.use = 'cc.pca',pt.size=2,group.by="cc.ident")
FeaturePlot(object = SO, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67","Cycs","Cdk1"), cols.use = c("blue", "yellow"), reduction.use = "umap",pt.size=1)
DimPlot(SO,reduction.use = 'umap', group.by="cc.ident",pt.size=2)

##### Atlas UMAP and markers visualization #####
genes.atlas <- rownames(SO.atlas@data)
SO.atlas <- FindVariableGenes(object = SO.atlas, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.05, x.high.cutoff = 4, y.cutoff = -0.1)
length(genes.atlas)
length(SO.atlas@var.genes)
SO.atlas@meta.data$group <- "atlas"
SO.atlas <- ScaleData(SO.atlas,do.scale=TRUE,do.center=TRUE,check.for.norm=TRUE)
SO.atlas <- RunPCA(SO.atlas,pcs.compute = 40)
PCElbowPlot(object = SO.atlas,num.pc = 40)
SO.atlas <- RunUMAP(SO.atlas,reduction.use="pca",dims.use=1:30,n_neighbors=10,min_dist=1,seed.use=655)
DimPlot(SO.atlas,reduction.use='umap',pt.size=1,do.label = TRUE)
tmp_markers <- c("Olfm4","Lgr5","Tff3","Fabp1","Mki67","Lyz1","Defa17","Chga","Krt19","Krt20","Dclk1","Anxa5")
tmp_markers <- c("nGene.log","nGene","percent.mito.log","percent.mito")
FeaturePlot(object = SO.atlas, features.plot = intersect(tmp_markers,genes.atlas), cols.use = c("blue", "yellow"), reduction.use = "umap",pt.size=1,no.axes = TRUE,no.legend=FALSE)

##### Atlas identities from the authors' metadata #####
atlas_names <- SO.atlas@cell.names
SO.atlas@ident <- factor(strsplit2(atlas_names,"_")[,3])
names(SO.atlas@ident) <- SO.atlas@cell.names
SO.atlas@meta.data$orig.cell.type <- SO.atlas@ident
DimPlot(SO.atlas,reduction.use = "umap")

##### Atlas cell cycle annotations #####
cc.genes <- readLines(con="../Resources/Regev_cellcyclegenes_mm.txt")
s.genes <- cc.genes[1:42]
g2m.genes <- cc.genes[43:93]
SO.atlas@var.genes <- union(SO.atlas@var.genes,intersect(cc.genes,rownames(SO.atlas@data)))
SO.atlas <- CellCycleScoring(object = SO.atlas, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
SO.atlas <- AddMetaData(object = SO.atlas, metadata = ordered(SO.atlas@ident, levels = c("G1","S","G2M")), col.name = "cc.ident")
# Plots on cell cycle: 
SO.atlas <- RunPCA(object = SO.atlas, pc.genes = c(s.genes, g2m.genes), do.print = FALSE,reduction.name="cc.pca")
DimPlot(object = SO.atlas,reduction.use = 'cc.pca',pt.size=2,group.by="cc.ident")
FeaturePlot(object = SO.atlas, features.plot = c("Pcna", "Top2a", "Mcm6", "Mki67","Cycs","Cdk1"), cols.use = c("blue", "yellow"), reduction.use = "umap",pt.size=1)
DimPlot(SO.atlas,reduction.use = 'umap', group.by="cc.ident",pt.size=2)

##### Atlas fuse identities in order to have consistent naming between the atlas and our datasets #####
SO.atlas@ident <- SO.atlas@meta.data$orig.cell.type
names(SO.atlas@ident) <- SO.atlas@cell.names
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Immature.Distal","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Immature.Proximal","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Mature.Distal","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Mature.Proximal","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Progenitor","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Progenitor.Early","G2/M")
SO.atlas <- RenameIdent(SO.atlas,"Enterocyte.Progenitor.Late","Enterocytes")
SO.atlas <- RenameIdent(SO.atlas,"Endocrine","Enteroendocrine")
SO.atlas <- RenameIdent(SO.atlas,"Stem","ISC")
SO.atlas <- RenameIdent(SO.atlas,"TA.Early","ISC")
SO.atlas <- RenameIdent(SO.atlas,"TA.G1","ISC")
SO.atlas <- RenameIdent(SO.atlas,"TA.G2","G2/M")
SO.atlas@ident %>% factor(levels = c("ISC","G2/M","Enterocytes","Paneth","Goblet","Tuft","Enteroendocrine")) -> SO.atlas@ident
SO.atlas@meta.data$cell.type <- SO.atlas@ident
DimPlot(SO.atlas,reduction.use = 'umap', group.by="ident", pt.size=2,cols.use = rainbow(length(levels(SO.atlas@ident))),do.label=TRUE)

