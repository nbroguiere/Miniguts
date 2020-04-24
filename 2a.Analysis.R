##### For exploration: Pick up one or two cluster to find DEGs between them / markers for user defined populations #####
ident.backup <- SO@ident
cluster.manual  <- DimPlot(SO,reduction.use = 'umap', do.identify=TRUE)
cluster.manual2 <- DimPlot(SO,reduction.use = 'umap', do.identify=TRUE)
ident.manual  <- ordered(ident.backup,c("Cluster1","Cluster2","Other"))
ident.manual[]<-"Other"
ident.manual[cluster.manual]  <- "Cluster1"
ident.manual[cluster.manual2] <- "Cluster2"
SO@ident <- ident.manual
DimPlot(SO,reduction.use = 'umap',pt.size=2)
markers <- FindMarkers(SO,ident.1="Cluster1",ident.2="Cluster2",only.pos = TRUE)
markers_sign <- markers[markers$p_val<0.05 & markers$avg_logFC>log(2),]

# # To save the ident:
# write.csv(x = SO@ident,file = "CustomPopulation.csv")
# SO@ident <- ident.backup

# # To restore a previous ident:
# tmp <- read.csv("./CustomPopulation.csv")
# SO@ident <- as.factor(tmp$x)
# names(SO@ident) <- names(ident_backup_groups)

writeClipboard(rownames(markers_sign)) # Paste in e.g. stringdb or pantherdb.org for functional annotations
writeClipboard(as.character(markers_sign$p_val)) # In case want to pick up the p-vals externally
writeClipboard(as.character(markers_sign$avg_logFC)) # Or the avglogFC

# Alternatively, directly find markers between two defined identity clusters and export as TSV:
markers <- FindMarkers(SO,ident.1="Enterocytes",ident.2="ISC",only.pos = TRUE)
write.table(markers, file="CustomMarkers.tsv", sep='\t',quote=TRUE, row.names=TRUE)

##### For exploration: Find and check features: ############
# Type it the root of a gene family of interest, and possibly edit the cells to use to plot only on one or another dataset:
Feat <- "Muc" 
tmp_markers <- grep(Feat,genes,value = T)
tmp_markers
FeaturePlot(object = SO, cells.use= c(cells.C1,cells.C2,cells.C3,cells.C4),features.plot = tmp_markers, cols.use = c("#EEEEEE", "blue"), reduction.use = "umap",pt.size = 1,no.legend = FALSE)

##### Cell counting by dataset and cell type #####
filename.tmp <- "Detailed_cell_counts.tsv"
cell.number.atlas <- matrix(nrow = length(levels(SO.atlas@meta.data$orig.ident)), ncol = length(levels(SO.atlas@meta.data$cell.type)))
rownames(cell.number.atlas) <- levels(SO.atlas@meta.data$orig.ident)
colnames(cell.number.atlas) <- levels(SO.atlas@meta.data$cell.type)
for(i in levels(SO.atlas@meta.data$orig.ident)){
  for(j in levels(SO.atlas@meta.data$cell.type)){
    cell.number.atlas[i,j] <- sum(SO.atlas@meta.data$orig.ident==i & SO.atlas@meta.data$cell.type==j)
  }
}
write.table(cell.number.atlas, file=filename.tmp, sep='\t',quote=FALSE, row.names=TRUE)

SO@meta.data$orig.ident <- factor(SO@meta.data$orig.ident)
cell.number <- matrix(nrow = length(levels(SO@meta.data$orig.ident)), ncol = length(levels(SO@meta.data$cell.type)))
rownames(cell.number) <- levels(SO@meta.data$orig.ident)
colnames(cell.number) <- levels(SO@meta.data$cell.type)
for(i in levels(SO@meta.data$orig.ident)){
  for(j in levels(SO@meta.data$cell.type)){
    cell.number[i,j] <- sum(SO@meta.data$orig.ident==i & SO@meta.data$cell.type==j)
  }
}
write.table(cell.number, file=filename.tmp, sep='\t',quote=FALSE, row.names=TRUE, append = TRUE)

##### Villus top enterocytes vs villus bot enterocytes: violin and scatter plots #####
ident_tmp <- as.character(ident_backup_groups)
ident_tmp[ident_tmp=="C3"] <- "Miniguts"
ident_tmp[ident_tmp=="C2"] <- "Miniguts"
ident_tmp[ident_tmp=="C1"] <- "Organoids"
SO@meta.data$simplified.ident <- factor(ident_tmp)
VlnPlot(object=SubsetData(SO,cells.use = c(cells.C1,cells.C2,cells.C3)),features.plot = "top.genes",ident.include = c("Enterocytes","TopEnterocytes"),group.by = "simplified.ident")

SO@ident <- factor(SO@meta.data$cell.type)
names(SO@ident) <- rownames(SO@meta.data)
SO.tmp <- SubsetData(SO,ident.use = c("Enterocytes","TopEnterocytes"))
SO.tmp@ident <- ident_backup_groups[colnames(SO.tmp@data)]
SO.tmp@ident[SO.tmp@ident=="C3"] <- "C2"
GenePlot(SO.tmp,"bot.genes","top.genes",cell.ids = intersect(colnames(SO.tmp@data),c(cells.C1,cells.C2,cells.C3)),cex.use = 0.5)

##### Defensins for Paneth cell highlight #####
#Check all defensins are markers of Paneth cells both in vitro and in vivo
defa.in.vitro <- grep("^Defa[0-9]*$",genes,value = T)
defa.in.atlas <- grep("^Defa[0-9]*$",genes.atlas,value = T)
defas.in.common <- sort(intersect(defa.in.vitro,defa.in.atlas))
FeaturePlot(object=SO.atlas,features.plot = defas.in.common,reduction.use = "umap", cols.use = c("lightgrey","blue"))
FeaturePlot(object=SO,features.plot = defas.in.common,reduction.use = "umap", cols.use = c("lightgrey","blue"))

# Defensins which are the most specific (least smear to Goblet) can be used to highlight Paneth cell position in UMAP:
defas.paneth  <- paste0("Defa",c(3,20,21,22,26))
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[defas.paneth,])), col.name = "defas.paneth")
SO.atlas <- AddMetaData(object = SO.atlas, metadata = colSums(as.matrix(SO.atlas@data[defas.paneth,])), col.name = "defas.paneth")
FeaturePlot(object=SO,features.plot = "defas.paneth",reduction.use = "umap", cols.use = c("lightgrey","blue"))
FeaturePlot(object=SO.atlas,features.plot = "defas.paneth",reduction.use = "umap", cols.use = c("lightgrey","blue"))

##### Violin and dot plots #####
# Matrix deposition
tmp_markers=c(paste0("Lama",1:5),paste0("Lamb",1:3),paste0("Lamc",1:3),"Nid1","Nid2","Col4a1","Col4a2","Hspg2","Fn1")
VlnPlot(object = SO, features.plot = intersect(tmp_markers,genes),y.max = 1.5,x.lab.rot = TRUE,size.x.use = 5, size.title.use = 9, point.size.use = 0.01,remove.legend = TRUE)

# Check gene correlations:
GenePlot(object = SO, "Lgr5", "Bmi1") # color by identity
ggplot(as.data.frame(t(as.matrix(SO@data)))) + geom_point(aes(x=Lama3, y=Lamb3, color=Lamc2))

# Dot plots of highlighted genes of interest for the microfold-like cells:
tmp_markers <- readLines(con="./Resources/Microfold_Highlights_v4.txt")
Mucins <- grep("Muc",genes,value=T)
SO <- AddMetaData(object = SO, metadata = colSums(as.matrix(SO@data[Mucins,])), col.name = "Mucins")
SO.tmp <- SO
SO.tmp@data <- rbind(SO.tmp@data,SO.tmp@meta.data$Mucins) # Create a virtual gene for mucins in a temporary Seurat object to seamlessly include this sum of genes in the dot plot. 
rownames(SO.tmp@data) <- c(rownames(SO@data),"Mucins")
DotPlot(object=SO.tmp,genes.plot = intersect(tmp_markers[length(tmp_markers):1],rownames(SO.tmp@data)), cols.use = c("lightgrey","red"),plot.legend = TRUE,col.min = 0,col.max = 1.5,x.lab.rot = TRUE) # Exported 25x3 inches
rm(SO.tmp)

##### Find cell type markers in vitro #####
# Compare each cluster to all the other clusters (cell types only, do not try to differ from G2/M), and then find the genes which are most significantly upregulated compared to the most similar cluster (i.e. max of p-values with other clusters)
markers <- NULL
for(i in levels(SO@ident)) # i is the cell type for which we look for markers
{
  for(j in setdiff(levels(SO@ident),c(i,"G2/M"))) # j is the cell type to which we compare. Do not compare any cells to G2/M. 
  {
    if((i=="ISC" & j=="Paneth")|(i=="Enterocytes" & j=="TopEnterocytes")|(i=="Goblet" & j=="Paneth")){} # Do not compare ISC to Paneth, enterocytes to Top enterocytes, because of natural expected shared markers. 
    else
    {
      markers[[i]][[j]] <- FindMarkers(object = SO, ident.1 = i, ident.2 = j, print.bar = TRUE,only.pos = TRUE,logfc.threshold = 0.2)
    }
  }
}
# Prefiltering for p value and log ratio:
markers_sign_tmp <- NULL
for(i in names(markers))
{
  for(j in names(markers[[i]]))
  {
    markers_sign_tmp[[i]][[j]] <- markers[[i]][[j]][markers[[i]][[j]]$p_val<0.05,]
    markers_sign_tmp[[i]][[j]] <- markers_sign_tmp[[i]][[j]][markers_sign_tmp[[i]][[j]]$avg_logFC>log(1.25),]
  }
}
# For further filtering (for heatmap only, keep all genes for significant markers with genes_tmp <- genes)
genes_tmp <- genes[-grep("^Rp",genes)] # Exclude ribosomal proteins
gene_dropout_rate <- rowSums(as.matrix(SO@data[genes,]==0))/length(colnames(SO@data))
genes_tmp <- intersect(genes_tmp,genes[gene_dropout_rate>0.25]) # Exclude genes which are not dropped out in at least 25% of the cells in the dataset

# Compute the minimal log FC compared to all other clusters, and maximal p-value. Currently use LogFC to pick the top genes for heatmap. 
markers_sign <- NULL
for(i in names(markers))
{
  tmp2 <- genes_tmp
  for(j in names(markers[[i]]))
  {
    tmp2 <- intersect(tmp2,rownames(markers_sign_tmp[[i]][[j]]))
  }
  markers_sign[[i]] <- data.frame(gene=tmp2,minlogFC=numeric(length(tmp2)),max_pval=numeric(length(tmp2)),stringsAsFactors=FALSE)
  rownames(markers_sign[[i]]) <- tmp2
  for(g in markers_sign[[i]]$gene)
  {
    avglogFC_tmp <- NULL
    p_val_tmp <- NULL
    for(j in names(markers[[i]]))
    {
      avglogFC_tmp <- c(avglogFC_tmp,markers_sign_tmp[[i]][[j]][g,]$avg_logFC)
      p_val_tmp <- c(p_val_tmp,markers_sign_tmp[[i]][[j]][g,]$p_val)
    }
    markers_sign[[i]][g,"minlogFC"]<- min(avglogFC_tmp)
    markers_sign[[i]][g,"max_pval"]<- max(p_val_tmp)
  }
}

top_conserved_markers <- NULL
n_genes_per_marker <- 20
for(i in levels(SO@ident))
{
  order_i <- order(-markers_sign[[i]]$minlogFC)
  markers_sign[[i]][order_i[1:n_genes_per_marker],] %>% rownames -> top_conserved_markers[[i]]
}

# Organize the cells by UMAP2 coordinate:
cells.use.invitro <- rownames(SO@dr$umap@cell.embeddings[order(-SO@dr$umap@cell.embeddings[,2]),])
DoHeatmap(object = SO, use.scaled = FALSE, genes.use = unlist(top_conserved_markers[]), slim.col.label = TRUE, remove.key = TRUE,col.low="red",col.mid="blue",col.high="yellow",draw.line = FALSE,title = NULL, cex.col = 5,cex.row = 6, group.label.loc = "top", group.label.rot = TRUE,group.cex = 10, assay.type = "RNA",do.plot = TRUE,cells.use = cells.use.invitro, disp.min = 0, disp.max = 5)

# Write all significant markers into a file:
write.table(x="CellTypeMarkersInVitro",file="CellTypeMarkersInVitro.tsv",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
for(i in names(markers_sign))
{
  write.table(x=i,file="CellTypeMarkersInVitro.tsv",append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(x=t(markers_sign[[i]]),file="CellTypeMarkersInVitro.tsv",append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

##### Find cell type markers in vivo  #####
# Compare each cluster to all the other clusters (cell types only, do not try to differ from G2/M), and then find the genes which are most significantly upregulated compared to the most similar cluster (i.e. max of p-values with other clusters)
markers <- NULL
for(i in levels(SO.atlas@ident))
{
  for(j in setdiff(levels(SO.atlas@ident),c(i,"G2/M"))) # Do not compare any cells to G2/M. 
  {
    if((i=="ISC" & j=="Paneth")|(i=="Goblet" & j=="Paneth")){} # Do not compare ISC to Paneth, because of natural expected shared markers.Do not compare Goblet to Paneth 
    else
    {
      markers[[i]][[j]] <- FindMarkers(object = SO.atlas, ident.1 = i, ident.2 = j, print.bar = TRUE,only.pos = TRUE,logfc.threshold = 0.2)
    }
  }
}
# Prefiltering for p value and log ratio:
markers_sign_tmp <- NULL
for(i in names(markers))
{
  for(j in names(markers[[i]]))
  {
    markers_sign_tmp[[i]][[j]] <- markers[[i]][[j]][markers[[i]][[j]]$p_val<0.05,]
    markers_sign_tmp[[i]][[j]] <- markers_sign_tmp[[i]][[j]][markers_sign_tmp[[i]][[j]]$avg_logFC>log(1.25),]
  }
}
# For further filtering (for heatmap only, keep all genes for significant markers with genes_tmp <- genes)
genes_tmp <- genes.atlas[-grep("^Rp",genes.atlas)] # Exclude ribosomal proteins
gene_dropout_rate <- rowSums(as.matrix(SO.atlas@data[genes.atlas,]==0))/length(colnames(SO.atlas@data))
genes_tmp <- intersect(genes_tmp,genes.atlas[gene_dropout_rate>0.25]) # Exclude genes which are not dropped out in at least 25% of the cells in the dataset

# Compute the minimal log FC compared to all other clusters, and maximal p-value. Use LogFC to select the top genes for heatmap. 
markers_sign <- NULL
for(i in names(markers))
{
  tmp2 <- genes_tmp
  for(j in names(markers[[i]]))
  {
    tmp2 <- intersect(tmp2,rownames(markers_sign_tmp[[i]][[j]]))
  }
  markers_sign[[i]] <- data.frame(gene=tmp2,minlogFC=numeric(length(tmp2)),max_pval=numeric(length(tmp2)),stringsAsFactors=FALSE)
  rownames(markers_sign[[i]]) <- tmp2
  for(g in markers_sign[[i]]$gene)
  {
    avglogFC_tmp <- NULL
    p_val_tmp <- NULL
    for(j in names(markers[[i]]))
    {
      avglogFC_tmp <- c(avglogFC_tmp,markers_sign_tmp[[i]][[j]][g,]$avg_logFC)
      p_val_tmp <- c(p_val_tmp,markers_sign_tmp[[i]][[j]][g,]$p_val)
    }
    markers_sign[[i]][g,"minlogFC"]<- min(avglogFC_tmp)
    markers_sign[[i]][g,"max_pval"]<- max(p_val_tmp)
  }
}

top_conserved_markers <- NULL
n_genes_per_marker <- 10
for(i in levels(SO.atlas@ident))
{
  order_i <- order(-markers_sign[[i]]$minlogFC)
  markers_sign[[i]][order_i[1:n_genes_per_marker],] %>% rownames -> top_conserved_markers[[i]]
}
cell.type.markers.atlas <- unlist(top_conserved_markers[])
DoHeatmap(object = SO.atlas, use.scaled = FALSE, genes.use = cell.type.markers.atlas, slim.col.label = TRUE, remove.key = TRUE,col.low="red",col.mid="blue",col.high="yellow",draw.line = FALSE,title = NULL, cex.col = 5,cex.row = 6, group.label.loc = "top", group.label.rot = TRUE,group.cex = 10, assay.type = "RNA",do.plot = TRUE,cells.use = cells.use.atlas, disp.min = 0, disp.max = 5, group.spacing = 0) # Export 20x14

# Write all significant markers into a file:
write.table(x="CellTypeMarkersInVivo",file="CellTypeMarkersInVivo.tsv",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
for(i in names(markers_sign))
{
  write.table(x=i,file="CellTypeMarkersInVivo.tsv",append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(x=t(markers_sign[[i]]),file="CellTypeMarkersInVivo.tsv",append=TRUE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
}

##### Heatmap in vitro from in vivo markers #####
DoHeatmap(object = SO, use.scaled = FALSE, genes.use = cell.type.markers.atlas, slim.col.label = TRUE, remove.key = TRUE,col.low="red",col.mid="blue",col.high="yellow",draw.line = FALSE,title = NULL, cex.col = 5,cex.row = 6, group.label.loc = "top", group.label.rot = TRUE,group.cex = 10, assay.type = "RNA",do.plot = TRUE,cells.use = cells.use.invitro, disp.min = 0, disp.max = 5, group.spacing = 0) # Export 20x14

##### Heatmap in vivo from in vitro markers #####
# Organize the cells by UMAP coordinate, to follow the stem-to-enterocyte axis:
cells.tmp <- colnames(SO.atlas@data)
cells.1 <- cells.tmp[SO.atlas@ident %in% c("ISC","Tuft","Endocrine","Goblet","Paneth")]
cells.12 <- cells.tmp[SO.atlas@ident=="G2/M"]
cells.2 <- cells.tmp[SO.atlas@ident=="Enterocytes"]

tmp <- SO.atlas@dr$umap@cell.embeddings[cells.1,]
cells.1  <- rownames(tmp[order(tmp[,1]),])
tmp <- SO.atlas@dr$umap@cell.embeddings[cells.12,]
cells.12  <- rownames(tmp[order(tmp[,1]-tmp[,2]),])
tmp <- SO.atlas@dr$umap@cell.embeddings[cells.2,]
cells.2  <- rownames(tmp[order(-tmp[,2]),])

cells.use.atlas <- c(cells.1,cells.12,cells.2)

# Load the gene list from the in vitro work: 
genes.use <- readLines("./Resources/HeatmapMarkerGeneListfromInVitro_noM.txt")
# Make the heatmap:
DoHeatmap(object = SO.atlas, use.scaled = FALSE, genes.use = genes.use, slim.col.label = TRUE, remove.key = TRUE,col.low="red",col.mid="blue",col.high="yellow",draw.line = FALSE,title = NULL, cex.col = 5,cex.row = 6, group.label.loc = "top", group.label.rot = TRUE,group.cex = 10, assay.type = "RNA",do.plot = TRUE,cells.use = cells.use.atlas, disp.min = 0, disp.max = 5, group.spacing = 10)

##### Isolate cell types and prepare data for GSEA in Java standalone software #####
# Documentation on:
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

enterocytes <- SubsetData(SO, ident.use = c("Enterocytes","TopEnterocytes"), subset.raw = T)
enterocytes <- SetAllIdent(enterocytes, id = "group")
ISC <- SubsetData(SO, ident.use = "ISC", subset.raw = T)
ISC <- SetAllIdent(ISC, id = "group")
Enteroendocrine <- SubsetData(SO, ident.use = "Enteroendocrine", subset.raw = T)
Enteroendocrine <- SetAllIdent(Enteroendocrine, id = "group")
Goblet <- SubsetData(SO, ident.use = "Goblet", subset.raw = T)
Goblet <- SetAllIdent(Goblet, id = "group")
Tuft <- SubsetData(SO, ident.use = "Tuft", subset.raw = T)
Tuft <- SetAllIdent(Tuft, id = "group")
Paneth <- SubsetData(SO, ident.use = "Paneth", subset.raw = T)
Paneth <- SetAllIdent(Paneth, id = "group")
Microfold <- SubsetData(SO, ident.use = "Microfold", subset.raw = T)
Microfold <- SetAllIdent(Microfold, id = "group")

celltype <- c("enterocytes","ISC","Enteroendocrine","Goblet","Tuft","Paneth","Microfold","AllCells")
target_folder <- "./InterferonResponse/"
for(i in 1:length(celltype))
{
  if(celltype[i]=="AllCells")
  {
    tmp <- SO
    filename <- celltype[i]
  } else {
    filename <- celltype[i]
    tmp <- eval(parse(text = celltype[i]))
  }
  cells_tmp <- intersect(c(cells.C2,cells.C4),colnames(tmp@data))
  tmp <- as.data.frame(as.matrix(tmp@data[,cells_tmp]))
  tmp2 <- matrix(data="na",nrow=nrow(tmp),ncol=1)
  tmp3 <- as.matrix(rownames(tmp))
  tmp <- cbind(tmp3,tmp2,tmp)
  #tmp <- rbind(colnames(tmp),tmp)
  #for .cgt files preface with: write.table(x=matrix(data=c("#1.2","",nrow(tmp),ncol(tmp)-2),nrow = 2,ncol=2,byrow = TRUE),file = paste(target_folder,filename,".txt",sep=""),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  write.table(x=t(as.matrix(c("NAME","Description",colnames(tmp)[3:length(colnames(tmp))]))),file = paste(target_folder,filename,".txt",sep=""),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  write.table(x=tmp,file = paste(target_folder,filename,".txt",sep=""),append = TRUE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  n_samples <- ncol(tmp)-2
  n_classes <- 2
  write.table(x=matrix(data=c(n_samples,n_classes,1,"#","C2","C4"),nrow = 2,ncol=3,byrow = TRUE),file = paste(target_folder,filename,".cls",sep=""),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  group_identity <- vector(mode="integer",length=length(cells_tmp))
  names(group_identity) <- cells_tmp
  group_identity[intersect(cells.C4,names(group_identity))] <- 1
  write.table(x=t(as.matrix(group_identity)),file = paste(target_folder,filename,".cls",sep=""),append = TRUE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
}

##### Plot genes up/downregulated upon treatment (volcano plots) #####
# Function definition
if(TRUE)
{
  LabelPoint <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                         adj.y.s = 0, text.size = 5, segment.size = 0.4) {
    for (i in genes) {
      x1 <- exp.mat[i, 1]
      y1 <- exp.mat[i, 2]
      plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                              label = i, size = text.size)
      plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                                adj.y.s, yend = y1, size = segment.size)
    }
    return(plot)
  }
  LabelPointRed <- function(plot, genes, exp.mat, adj.x.t = 0, adj.y.t = 0, adj.x.s = 0, 
                            adj.y.s = 0, text.size = 5, segment.size = 0.4) {
    for (i in genes) {
      x1 <- exp.mat[i, 1]
      y1 <- exp.mat[i, 2]
      plot <- plot + annotate("text", x = x1 + adj.x.t, y = y1 + adj.y.t, 
                              label = i, size = text.size, colour = "red")
      plot <- plot + annotate("segment", x = x1 + adj.x.s, xend = x1, y = y1 + 
                                adj.y.s, yend = y1, size = segment.size)
    }
    return(plot)
  }
  
  LabelUR <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.r.t = 0.15, adj.u.s = 0.05, 
                      adj.r.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = adj.r.t, 
                      adj.y.s = adj.u.s, adj.x.s = adj.r.s, ...))
  }
  
  LabelUL <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                      adj.l.s = 0.05, ...) {
    return(LabelPoint(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                      adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
  }
  LabelULRed <- function(plot, genes, exp.mat, adj.u.t = 0.1, adj.l.t = 0.15, adj.u.s = 0.05, 
                         adj.l.s = 0.05, ...) {
    return(LabelPointRed(plot, genes, exp.mat, adj.y.t = adj.u.t, adj.x.t = -adj.l.t, 
                         adj.y.s = adj.u.s, adj.x.s = -adj.l.s, ...))
  }
}

# Subset the cell types and compute the differential expression:
CA <- "C2" # First condition to compare
CB <- "C4" # Second condition to compare

# Use both custom highlighted interferon related genes and the MSigDB interferon alpha hallmark response converted to mm based this db for conversion: http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
Ifa_resp <- read.table(file="./Resources/mm_GSEA_HALLMARK_INTERFERON_ALPHA_RESPONSE.txt",header=FALSE,quote="")
Ifa_resp <- data.frame(lapply(Ifa_resp, as.character), stringsAsFactors=FALSE)$V1
custom_Ifa_genes <- c("Ifitm3","Ifi27","Ifi27l2b","Bst2","B2m","Lgals9","Tcim","Ifnar1","Ddx58","Oasl1","Oasl2","Ifit3","Ifit3b","Ifit1","Ifit2","Ifit3","Ifitm1","Ifitm2","Ifitm3","Oasl1","Oasl2")
Ifa_resp <- c(Ifa_resp,custom_Ifa_genes)

# For plotting volcano plots (-log(Pvalue) vs log ratio)
p_val_thresh <- 0.05
LR_thresh <- log(1.25) # Used 1.5 for paper

enterocytes_markers <- FindMarkers(enterocytes,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
enterocytes_markers$p_val_mlog10 <- -log10(enterocytes_markers$p_val)
enterocytes_markers$genes <- rownames(enterocytes_markers)
enterocytes_markers_significant <- enterocytes_markers[enterocytes_markers$p_val<p_val_thresh & enterocytes_markers$avg_logFC>LR_thresh,]
enterocytes_markers_significant_down <- enterocytes_markers[enterocytes_markers$p_val<p_val_thresh & enterocytes_markers$avg_logFC<(-LR_thresh),];

ISC_markers <- FindMarkers(ISC,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
ISC_markers$p_val_mlog10 <- -log10(ISC_markers$p_val)
ISC_markers$genes <- rownames(ISC_markers)
ISC_markers_significant <- ISC_markers[ISC_markers$p_val<p_val_thresh & ISC_markers$avg_logFC>LR_thresh,]
ISC_markers_significant_down <- ISC_markers[ISC_markers$p_val<p_val_thresh & ISC_markers$avg_logFC<(-LR_thresh),]

Enteroendocrine_markers <- FindMarkers(Enteroendocrine,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
Enteroendocrine_markers$p_val_mlog10 <- -log10(Enteroendocrine_markers$p_val)
Enteroendocrine_markers$genes <- rownames(Enteroendocrine_markers)
Enteroendocrine_markers_significant <- Enteroendocrine_markers[Enteroendocrine_markers$p_val<p_val_thresh & Enteroendocrine_markers$avg_logFC>LR_thresh,]
Enteroendocrine_markers_significant_down <- Enteroendocrine_markers[Enteroendocrine_markers$p_val<p_val_thresh & Enteroendocrine_markers$avg_logFC<(-LR_thresh),]

Goblet_markers <- FindMarkers(Goblet,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
Goblet_markers$p_val_mlog10 <- -log10(Goblet_markers$p_val)
Goblet_markers$genes <- rownames(Goblet_markers)
Goblet_markers_significant <- Goblet_markers[Goblet_markers$p_val<p_val_thresh & Goblet_markers$avg_logFC>LR_thresh,]
Goblet_markers_significant_down <- Goblet_markers[Goblet_markers$p_val<p_val_thresh & Goblet_markers$avg_logFC<(-LR_thresh),]

Tuft_markers <- FindMarkers(Tuft,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
Tuft_markers$p_val_mlog10 <- -log10(Tuft_markers$p_val)
Tuft_markers$genes <- rownames(Tuft_markers)
Tuft_markers_significant <- Tuft_markers[Tuft_markers$p_val<p_val_thresh & Tuft_markers$avg_logFC>LR_thresh,]
Tuft_markers_significant_down <- Tuft_markers[Tuft_markers$p_val<p_val_thresh & Tuft_markers$avg_logFC<(-LR_thresh),]

Paneth_markers <- FindMarkers(Paneth,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
Paneth_markers$p_val_mlog10 <- -log10(Paneth_markers$p_val)
Paneth_markers$genes <- rownames(Paneth_markers)
Paneth_markers_significant <- Paneth_markers[Paneth_markers$p_val<p_val_thresh & Paneth_markers$avg_logFC>LR_thresh,]
Paneth_markers_significant_down <- Paneth_markers[Paneth_markers$p_val<p_val_thresh & Paneth_markers$avg_logFC<(-LR_thresh),]

Microfold_markers <- FindMarkers(Microfold,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
Microfold_markers$p_val_mlog10 <- -log10(Microfold_markers$p_val)
Microfold_markers$genes <- rownames(Microfold_markers)
Microfold_markers_significant <- Microfold_markers[Microfold_markers$p_val<p_val_thresh & Microfold_markers$avg_logFC>LR_thresh,]
Microfold_markers_significant_down <- Microfold_markers[Microfold_markers$p_val<p_val_thresh & Microfold_markers$avg_logFC<(-LR_thresh),]

p1 <- ggplot(enterocytes_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Enterocytes")
p1 <- p1 + geom_point(colour="red",alpha=1/10,size=1)
p1 <- p1 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(enterocytes_markers_significant))
p1 <- LabelULRed(p1, genes = gtmp, enterocytes_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(enterocytes_markers_significant$genes,gtmp)
#p1 <- LabelUL(p1, genes = gtmp, enterocytes_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p2 <- ggplot(ISC_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("ISC and TA")
p2 <- p2 + geom_point(colour="red",alpha=1/10,size=1)
p2 <- p2 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(ISC_markers_significant))
p2 <- LabelULRed(p2, genes = gtmp, ISC_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(ISC_markers_significant$genes,gtmp)
#p2 <- LabelUL(p2, genes = gtmp, ISC_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p3 <- ggplot(Enteroendocrine_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Enteroendocrine")
p3 <- p3 + geom_point(colour="red",alpha=1/10,size=1)
p3 <- p3 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(Enteroendocrine_markers_significant))
p3 <- LabelULRed(p3, genes = gtmp, Enteroendocrine_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(Enteroendocrine_markers_significant$genes,gtmp)
#p3 <- LabelUL(p3, genes = gtmp, Enteroendocrine_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p4 <- ggplot(Goblet_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Goblet cells")
p4 <- p4 + geom_point(colour="red",alpha=1/10,size=1)
p4 <- p4 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(Goblet_markers_significant))
p4 <- LabelULRed(p4, genes = gtmp, Goblet_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(Goblet_markers_significant$genes,gtmp)
#p4 <- LabelUL(p4, genes = gtmp, Goblet_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p5 <- ggplot(Tuft_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Tuft cells")
p5 <- p5 + geom_point(colour="red",alpha=1/10,size=1)
p5 <- p5 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(Tuft_markers_significant))
p5 <- LabelULRed(p5, genes = gtmp, Tuft_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(Tuft_markers_significant$genes,gtmp)
#p5 <- LabelUL(p5, genes = gtmp, Tuft_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p6 <- ggplot(Paneth_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Paneth cells")
p6 <- p6 + geom_point(colour="red",alpha=1/10,size=1)
p6 <- p6 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(Paneth_markers_significant))
p6 <- LabelULRed(p6, genes = gtmp, Paneth_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(Paneth_markers_significant$genes,gtmp)
#p6 <- LabelUL(p6, genes = gtmp, Paneth_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

p7 <- ggplot(Microfold_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("Microfold cells")
p7 <- p7 + geom_point(colour="red",alpha=1/10,size=1)
p7 <- p7 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(Microfold_markers_significant))
p7 <- LabelULRed(p7, genes = gtmp, Microfold_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)
#gtmp <- setdiff(Microfold_markers_significant$genes,gtmp)
#p7 <- LabelUL(p7, genes = gtmp, Microfold_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 0.2,adj.l.t = -0.2,adj.u.s = 0.075,adj.l.s = -0.075)

plot_grid(p1,p2,p3,p4,p5,p6)

write.table(enterocytes_markers_significant, file="enterocytes_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)
write.table(ISC_markers_significant, file="ISC_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)
write.table(Enteroendocrine_markers_significant, file="Enteroendocrine_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)
write.table(Goblet_markers_significant, file="Goblet_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)
write.table(Tuft_markers_significant, file="Tuft_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)
write.table(Paneth_markers_significant, file="Paneth_markers_significant.tsv", sep='\t',quote=TRUE, row.names=TRUE)

writeClipboard(as.character(Paneth_markers$avg_logFC)) # For copy-pasting into excel sheet or databases like interferome and string-db for functional annotations

# Similar analysis with all the cells pulled together:
CA <- "C2" # First condition to compare
CB <- "C4" # Second condition to compare

SO_tmp <- SO
SO_tmp <- SetAllIdent(SO_tmp,id="group")

all_cells_markers <- FindMarkers(SO_tmp,ident.1=CB,ident.2=CA,logfc.threshold=0,test.use="wilcox",min.pct=0,min.cells.gene=1,min.cells.group=1)
all_cells_markers$p_val_mlog10 <- -log10(all_cells_markers$p_val)
all_cells_markers$genes <- rownames(all_cells_markers)
p_val_thresh <- 0.01
LR_thresh <- log(1.25)
all_cells_markers_significant <- all_cells_markers[all_cells_markers$p_val<p_val_thresh & all_cells_markers$avg_logFC>LR_thresh,]
all_cells_markers_significant_down <- all_cells_markers[all_cells_markers$p_val<p_val_thresh & all_cells_markers$avg_logFC<(-LR_thresh),]
rm(SO_tmp)

p7 <- ggplot(all_cells_markers,aes_string("avg_logFC", "p_val_mlog10")) + geom_point(colour="blue",alpha=1,size=1) + ggtitle("All cells")
p7 <- p7 + geom_point(colour="red",alpha=1/10,size=1)
p7 <- p7 + geom_point(colour="yellow",alpha=1/100,size=1)
gtmp <- intersect(Ifa_resp,rownames(all_cells_markers))
p7 <- LabelULRed(p7, genes = gtmp, all_cells_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 4,adj.l.t = -0.1,adj.u.s = 4,adj.l.s = -0.1)
gtmp <- setdiff(all_cells_markers_significant$genes,gtmp)
p7 <- LabelUL(p7, genes = gtmp, all_cells_markers[gtmp,c("avg_logFC", "p_val_mlog10","genes")],adj.u.t = 4,adj.l.t = -0.1,adj.u.s = 4,adj.l.s = -0.1)

plot_grid(p7)
writeClipboard(as.character(all_cells_markers$genes)) # For pasting in string-db for functional annotations

#backup_C2vsC4 <- all_cells_markers
#backup_C2vsC4_significant <- all_cells_markers_significant

significant_genes_filtered <- setdiff(all_cells_markers_significant$genes,backup_C2vsC4_significant$genes)
writeClipboard(as.character(backup_C2vsC4$avg_logFC)) # For pasting in string-db for functional annotations

##### Find markers between treatments #####
Compare_1 <- "C2 Paneth"
Compare_2 <- "C4 Paneth"
ident.backup <- SO@ident
ident.manual <- ordered(paste(ident_backup_groups,ident_backup_clusters))
names(ident.manual) <- names(ident.backup)
SO@ident <- ident.manual
DimPlot(SO,reduction.use = 'umap',pt.size=2)
markers <- FindMarkers(SO,ident.1=Compare_1,ident.2=Compare_2,only.pos = FALSE,test.use="tobit")
writeClipboard(rownames(markers)) # Paste in pantherdb.org for functional annotations
SO@ident <- ident.backup

write.table(markers, file="C2vsC4_Paneth_Markers_tobit.tsv", sep='\t',quote=TRUE, row.names=TRUE)

