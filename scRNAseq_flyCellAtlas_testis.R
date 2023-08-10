# Script name:    scRNAseq_flyCellAtlas_testis.R
# Created on:     August_10_2023
# Author:         Dr. Martin Calvino
# Purpose:        Analyse cell-type expression patterns of NXF family members in Drosophila melanogaster's testis
        
# load libraries
library(SeuratDisk)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(scales)
library(MetBrewer)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork)

# connect to loom file
my.loom.file = file.choose()
my.loom.file # "/Users/martincalvinotorterolo/Desktop/myR_scripts/data/r_fca_biohub_testis_10x.loom"
testis.loom <- Connect(filename = my.loom.file, mode = "r+")
testis.loom

# interact with loom object and view matrix dataset
testis.loom[["attrs"]]
testis.loom[["col_attrs"]]
testis.loom[["row_attrs"]]
testis.loom[["matrix"]]

# Access the upper left corner of the data matrix
testis.loom[["matrix"]][1:15, 1:15]

# list the attributes of testis.loom object
attributes(testis.loom) # "loom"       "scdisk"     "H5File"     "H5RefClass" "R6"
class(testis.loom) # "loom"       "scdisk"     "H5File"     "H5RefClass" "R6"

# now convert H5AD file to Seurat (since .loom to Seurat does not work)
my.h5ad.file = file.choose()
my.h5ad.file # "/Users/martincalvinotorterolo/Desktop/myR_scripts/data/r_fca_biohub_testis_10x.h5ad"
Convert(my.h5ad.file, dest = "h5seurat", overwrite = TRUE)
seurat_testis <- LoadH5Seurat("/Users/martincalvinotorterolo/Desktop/myR_scripts/data/r_fca_biohub_testis_10x.h5seurat")
seurat_testis # seurat object from h5ad file contains only 2141 features as opposed to 15,833 features in the loom file

seurat_testis@misc$ID_categories #suerat_testis@misc is where the relevant info is

# reconstruct seurat object from loom and h5ad files
fca_testis_mat <- testis.loom[["/matrix"]][,]
fca_testis_mat <- Matrix::Matrix(fca_testis_mat, sparse = TRUE)
fca_testis_mat <- Matrix::t(fca_testis_mat)
fca_testis_mat # now I have genes as rows and cells as columns

# now extract the cell id and gene id and set them as the column and row names of the matrix
fca_testis_cellid <- testis.loom[["/col_attrs/CellID"]][]
fca_testis_cellid # example > "AAACCCACATTACTCT-6e669170__FCA59_Male_testis_adult_1dWT_Fuller_sample1"
fca_testis_geneid <- testis.loom[["row_attrs/Gene"]][]
fca_testis_geneid # example > "Atg1"

colnames(fca_testis_mat) <- fca_testis_cellid
rownames(fca_testis_mat) <- fca_testis_geneid
fca_testis_mat

# Pull three bits of metadata from the column attributes
attrs <- c("CellID", "ClusterID", "n_counts", "n_genes", "percent_mito", "annotation", "annotation_broad", "sex")
attrs.df <- map_dfc(attrs, ~ testis.loom[[paste0("col_attrs/", .)]][]) %>% as.data.frame()
colnames(attrs.df) <- attrs
rownames(attrs.df) <- fca_testis_cellid
View(attrs.df)

# build seurat object
fca_testis_seurat <- CreateSeuratObject(counts = fca_testis_mat, meta.data = attrs.df)
fca_testis_seurat
head(fca_testis_seurat, n=25)

fca_testis_h5 <- LoadH5Seurat("/Users/martincalvinotorterolo/Desktop/myR_scripts/data/r_fca_biohub_testis_10x.h5seurat")
fca_testis_seurat@reductions <- fca_testis_h5@reductions
rownames(fca_testis_seurat@reductions$pca@cell.embeddings) <- fca_testis_cellid
rownames(fca_testis_seurat@reductions$tsne@cell.embeddings) <- fca_testis_cellid
rownames(fca_testis_seurat@reductions$umap@cell.embeddings) <- fca_testis_cellid

fca_testis_seurat

# Visualize UMAP with clusters assigned to cell types and expression of NXF family members
DimPlot(fca_testis_seurat, reduction = "umap", group.by = 'annotation', label = T, label.size = 2.5) + NoLegend()
FeaturePlot(fca_testis_seurat, reduction = 'umap', features = c("sbr", "nxf2", "Nxf3", "nxf4"))

# Visualize percentage of cells in each cell types expressing NXF family members
DotPlot(fca_testis_seurat, features = c("sbr", "nxf2", "Nxf3", "nxf4", "mael",
                                        "piwi", "del", "BoYb", "shu", "zuc", "aub", "tej", "krimp", "vls", "Nbr", "csul",
                                        "cuff", "qin", "egg", "vret", "Nxt1", "arx", "Hen1", "Panx", "mino", "spn-E",
                                        "AGO3", "rhi", "Gasz", "Su(var)3-3", "wde", "SoYb", "armi"), group.by = 'annotation') + RotatedAxis()

NXF.testis <- DotPlot(fca_testis_seurat, features = c("sbr", "nxf2", "Nxf3", "nxf4", "mael","piwi", "del", 
                                                      "BoYb", "shu", "zuc", "aub", "tej", "krimp", "vls", 
                                                      "Nbr", "csul", "cuff", "qin", "egg", "vret", "Nxt1", 
                                                      "arx", "Hen1", "Panx", "mino", "spn-E", "AGO3", "rhi", 
                                                      "Gasz", "Su(var)3-3", "wde", "SoYb", "armi"), 
                      group.by = 'annotation') + RotatedAxis()

NXF.testis$data
Average_Expression_Scaled <- NXF.testis$data$avg.exp.scaled
Percent_Expressed <- NXF.testis$data$pct.exp

id <- factor(NXF.testis$data$id, order = TRUE, levels = c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition",
  "spermatocyte",
  "spermatocyte 0",
  "spermatocyte 1",
  "spermatocyte 2",
  "spermatocyte 3",
  "spermatocyte 4",
  "spermatocyte 5",
  "spermatocyte 6",
  "spermatocyte 7",
  "spermatocyte 7a",
  "late primary spermatocyte",
  "spermatid",
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid",
  "cyst stem cell",
  "early cyst cell 1",
  "early cyst cell 2",
  "spermatocyte cyst cell branch A",
  "spermatocyte cyst cell branch B",
  "cyst cell branch b"
))

ggplot(data = NXF.testis$data, aes(x = factor(NXF.testis$data$id, order = TRUE, levels = c(
  "spermatogonium",
  "mid-late proliferating spermatogonia",
  "spermatogonium-spermatocyte transition",
  "spermatocyte",
  "spermatocyte 0",
  "spermatocyte 1",
  "spermatocyte 2",
  "spermatocyte 3",
  "spermatocyte 4",
  "spermatocyte 5",
  "spermatocyte 6",
  "spermatocyte 7",
  "spermatocyte 7a",
  "late primary spermatocyte",
  "spermatid",
  "early elongation stage spermatid",
  "early-mid elongation-stage spermatid",
  "mid-late elongation-stage spermatid",
  "cyst stem cell",
  "early cyst cell 1",
  "early cyst cell 2",
  "spermatocyte cyst cell branch A",
  "spermatocyte cyst cell branch B",
  "cyst cell branch b")), 
  y = NXF.testis$data$features.plot,
  size = Percent_Expressed, fill = Average_Expression_Scaled)) +
  geom_point(shape = 21, colour = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust =1, vjust = 1)) +
  scale_fill_gradientn(colors = met.brewer('Derain')) +
  labs(x = "", y = "") +
  scale_size_area()

# all genes
markers <- NXF.testis$data$features.plot %>% 
  unique()

View(markers)

mat <- NXF.testis$data %>%
  filter(features.plot %in% markers) %>%
  select(-avg.exp, -pct.exp) %>%
  pivot_wider(names_from=id, values_from=avg.exp.scaled) %>%
  data.frame()

View(mat)

row.names(mat) <- mat$features.plot
View(mat)

mat <- mat[, -1]
View(mat)

clust <- hclust(dist(mat %>% as.matrix()))
plot(clust)

ddgram <- as.dendrogram(clust)
ggtree_plot <- ggtree::ggtree(ddgram)
ggtree_plot

my.dotPlot <- NXF.testis$data %>% 
  filter(features.plot %in% markers) %>%
  filter(pct.exp > 1) %>%
  mutate(features.plot = factor(features.plot, levels=clust$labels[clust$order])) %>%
  ggplot(aes(x=factor(id, order=TRUE, levels=c("germinal proliferation center hub", 
                                               "spermatogonium",
                                               "mid-late proliferating spermatogonia",
                                               "spermatogonium-spermatocyte transition",
                                               "spermatocyte",
                                               "spermatocyte 0",
                                               "spermatocyte 1",
                                               "spermatocyte 2",
                                               "spermatocyte 3",
                                               "spermatocyte 4",
                                               "spermatocyte 5",
                                               "spermatocyte 6",
                                               "spermatocyte 7",
                                               "spermatocyte 7a",
                                               "maturing primary spermatocyte",
                                               "late primary spermatocyte",
                                               "spermatid",
                                               "early elongation stage spermatid",
                                               "early-mid elongation-stage spermatid",
                                               "mid-late elongation-stage spermatid",
                                               "cyst stem cell",
                                               "early cyst cell 1",
                                               "early cyst cell 2",
                                               "spermatocyte cyst cell branch a",
                                               "spermatocyte cyst cell branch b",
                                               "cyst cell branch b",
                                               "cyst cell branch a",
                                               "cyst cell intermediate",
                                               "head cyst cell",
                                               "late cyst cell branch b",
                                               "late cyst cell branch a",
                                               "adult tracheocyte",
                                               "adult neuron",
                                               "adult fat body",
                                               "hemocyte",
                                               "muscle cell",
                                               "pigment cell",
                                               "secretory cell of the male reproductive tract",
                                               "male gonad associated epithelium",
                                               "unannotated"
  )), y=features.plot, fill=avg.exp.scaled, size=pct.exp)) +
  cowplot::theme_cowplot() +
  geom_point(shape=21, colour="black") + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) + 
  scale_fill_gradientn(colors=met.brewer('Derain')) + labs(x="", y="") + 
  scale_size_area()

plot_grid(ggtree_plot, NULL, my.dotPlot, nrow=1, rel_widths = c(0.25, -0.005, 2), align= 'h') 