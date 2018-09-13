#!/usr/bin/Rscript
library(getopt)
library(Seurat)
library(Matrix)

Func_OutUsage = function(x){
write("
Description:
	Identify clusters of cells by Seurat
", stdout());
write(getopt(x, usage=TRUE, command = paste("\nRscript ", get_Rscript_filename())) ,stdout());
write("
Note:
	'-t' option, a sample name per line, without title
Examples:
	Rscript ClusterBySeurat.R -i in_mat.txt -r rm.list -o out.prefix
", stdout())
q(status=1);
}
ParaSet = matrix( c(
	"input", "i", 1, "character", "* input matrix",
	"rmlist", "t", 2, "character", "sample names that should be removed",
	"dim", "n", 2, "character", "'dims.use' in 'Seurat::FindClusters', default 1:6",
	"k.param", "k", 2, "numerical", "'k.param' in 'Seurat::FindClusters'", 
	"resolution", "r", 2, "numerical", "'resolution' in 'Seurat::FindClusters'", 
	"out.prefix", "o", 2, "character", "output prefix"
	),byrow=TRUE, ncol=5
);
Para = getopt(ParaSet);
if (is.null(Para$input) || !(file.exists(Para$input))) { Func_OutUsage(ParaSet);}
if (is.null(Para$rmlist)) { Para$rmlist <- ''; }
if (is.null(Para$out.prefix)) { Para$out.prefix = "Output"; }


exp <- read.table(Para$input, header = T, check.names = F, stringsAsFactors = T)
if (Para$rmlist != ""){
 	rm <- read.table(Para$rmlist,header = F)
 	exp <- exp[apply(exp, 1, function(x) sum(x>1)>3), !colnames(exp)%in%rm$V1]
}else{
	exp <- exp[apply(exp, 1, function(x) sum(x>1)>3), ]
}

Obj <- CreateSeuratObject(raw.data = exp)
Obj <- NormalizeData(Obj)
Obj <- FindVariableGenes(object = Obj, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.5, x.high.cutoff = 10, y.cutoff = 0.5)

pdf(paste0(Para$out.prefix, ".pc_select.pdf"),width = 10,height = 10)
Obj <- ScaleData(object = Obj)
Obj <- RunPCA(object = Obj, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
Obj <- ProjectPCA(object = Obj, do.print = FALSE)
PCAPlot(object = Obj, dim.1 = 1, dim.2 = 2,pt.size = 1)
PCHeatmap(Obj, pc.use = 1:15, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
Obj <- JackStraw(Obj, num.replicate = 100)
JackStrawPlot(Obj, PCs = 1:20)
PCElbowPlot(Obj)
dev.off()

pdf(paste0(Para$out.prefix, ".tsne_cluster.pdf"),width = 8,height = 7)
Obj <- FindClusters(Obj, reduction.type = "pca",dims.use = c(1:6), save.SNN = T, 
                    k.param = 6, resolution = 0.6)

Obj <- RunTSNE(Obj, dims.use = c(1:6), do.fast = T, perplexity = 5)
TSNEPlot(Obj,pt.size = 3, do.label = F)
PCAPlot(object = Obj, dim.1 = 1, dim.2 = 2,pt.size = 1)
dev.off()

cluster.markers <- FindAllMarkers(Obj, only.pos = TRUE,  min.pct = 0.5, logfc.threshold = 0.1)
filter <- cluster.markers[cluster.markers$p_val<0.05, ]
pdf(paste0(Para$out.prefix, ".marker.heatmap.pdf"), 8, 12)
DoHeatmap(Obj, genes.use = filter$gene, slim.col.label = TRUE, remove.key = FALSE, cex.row=2)
dev.off()
write.table(filter,file = paste0(Para$out.prefix, ".cluster_markers.txt"), sep = "\t", quote = F)

write.table(Obj@dr$tsne@cell.embeddings, file = paste0(Para$out.prefix, ".tsne_coordinate.txt"), sep = "\t", quote = F)
cell_anno <- data.frame(Obj@ident)
write.table(cell_anno, file = paste0(Para$out.prefix, ".cluster.txt"), sep = "\t", quote = F)

save.image(paste0(Para$out.prefix, ".Seurat.RData"))
