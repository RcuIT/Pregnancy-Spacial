
library(Seurat)
library(scCustomize) 
library(ComplexHeatmap)

E_MTAB_6678_data <- read.delim("raw_data_ss2.txt",sep = "\t", header = T,check.names = FALSE)

rownames(E_MTAB_6678_data) <-E_MTAB_6678_data$Gene

E_MTAB_6678_data <- E_MTAB_6678_data[,-1]

E_MTAB_6678_obj <- CreateSeuratObject(E_MTAB_6678_data)
rownames(E_MTAB_6678_obj)
colnames(E_MTAB_6678_obj)

E_MTAB_6678_obj[[]]

# load the metada
E_MTAB_6678_metadata <- read.delim("meta_ss2.txt",sep = "\t", header = T,check.names = FALSE)

# check the metadata
head(E_MTAB_6678_metadata)

# Check for identical rownames between two objects
all.equal(rownames(E_MTAB_6678_obj[[]]), rownames(E_MTAB_6678_metadata))

# this was not true. the rows in the metadata and seurat metadata rows are not in the same order
# thus i have to rearrange the metadata rows as in the order of the seurat object

# re-arrangeging metadata rownames according the seurat obj rownames
E_MTAB_6678_metadata <- E_MTAB_6678_metadata[rownames(E_MTAB_6678_obj[[]]),
                                             , drop = FALSE]

# Check for identical rownames between two objects
all.equal(rownames(E_MTAB_6678_obj[[]]), rownames(E_MTAB_6678_metadata))

# now the rows in the both data set are TRUE , ie in the same order

# now transfer the metadata to seurat obj

E_MTAB_6678_obj <- AddMetaData(E_MTAB_6678_obj, 
                               metadata = E_MTAB_6678_metadata)
# check the metadata
E_MTAB_6678_obj[[]]

saveRDS(E_MTAB_6678_obj, "E_MTAB_6678_obj.rds")


########################

E_MTAB_6701_data <- read.delim("raw_data_10x.txt",sep = "\t", header = T,check.names = FALSE)
head(E_MTAB_6701_data)
rownames(E_MTAB_6701_data) <-E_MTAB_6701_data$Gene

E_MTAB_6701_data <- E_MTAB_6701_data[,-1]
E_MTAB_6701_data[1:2,1:2]

gene_symbols <- sub("_ENSG.*", "", rownames(E_MTAB_6701_data))
rownames(E_MTAB_6701_data) <- make.unique(gene_symbols)
E_MTAB_6701_data[1:2,1:2]

E_MTAB_6701_obj <- CreateSeuratObject(E_MTAB_6701_data)
rownames(E_MTAB_6701_obj)
colnames(E_MTAB_6701_obj)

length(E_MTAB_6701_obj@meta.data$orig.ident)
obj.metadata <- rownames(E_MTAB_6701_obj@meta.data)
length(obj.metadata)
head(E_MTAB_6701_obj@meta.data)

# load the metada
E_MTAB_6701_metadata <- read.delim("meta_10x.txt",sep = "\t", header = T,check.names = FALSE)
head(E_MTAB_6701_metadata)

metadata <- rownames(E_MTAB_6701_metadata)
length(metadata)

length(E_MTAB_6701_metadata$Fetus)
head(E_MTAB_6701_metadata)

# check the metadata
head(E_MTAB_6701_metadata)

# Check for identical rownames between two objects
all.equal(obj.metadata, metadata)

# this was not true. the rows in the metadata and seurat metadata rows are not in the same order
# thus i have to rearrange the metadata rows as in the order of the seurat object

# re-arrangeging metadata rownames according the seurat obj rownames
E_MTAB_6701_metadata <- E_MTAB_6701_metadata[rownames(E_MTAB_6701_obj[[]]),
                                             , drop = FALSE]

# Check for identical rownames between two objects
all.equal(rownames(E_MTAB_6701_obj[[]]), rownames(E_MTAB_6701_metadata))

# now the rows in the both data set are TRUE , ie in the same order

# now transfer the metadata to seurat obj

E_MTAB_6701_obj <- AddMetaData(E_MTAB_6701_obj, 
                               metadata = E_MTAB_6701_metadata)
# check the metadata
E_MTAB_6701_obj[[]]

saveRDS(E_MTAB_6701_obj, "E_MTAB_6701_obj2.rds")

####################################

E_MTAB_6701_data <- read.delim("raw_data_10x.txt",sep = "\t", header = T,check.names = FALSE)

rownames(E_MTAB_6701_data) <-E_MTAB_6701_data$Gene

E_MTAB_6701_data <- E_MTAB_6701_data[,-1]

E_MTAB_6701_obj <- CreateSeuratObject(E_MTAB_6701_data)
rownames(E_MTAB_6701_obj)
colnames(E_MTAB_6701_obj)

E_MTAB_6701_obj[[]]

# load the metada
E_MTAB_6701_metadata <- read.delim("meta_10x.txt",sep = "\t", header = T,check.names = FALSE)

# check the metadata
head(E_MTAB_6701_metadata)

# Check for identical rownames between two objects
all.equal(rownames(E_MTAB_6701_obj[[]]), rownames(E_MTAB_6701_metadata))

# this was not true. the rows in the metadata and seurat metadata rows are not in the same order
# thus i have to rearrange the metadata rows as in the order of the seurat object

# re-arrangeging metadata rownames according the seurat obj rownames
E_MTAB_6701_metadata <- E_MTAB_6701_metadata[rownames(E_MTAB_6701_obj[[]]),
                                             , drop = FALSE]

# Check for identical rownames between two objects
all.equal(rownames(E_MTAB_6701_obj[[]]), rownames(E_MTAB_6701_metadata))

# now the rows in the both data set are TRUE , ie in the same order

# now transfer the metadata to seurat obj

E_MTAB_6701_obj <- AddMetaData(E_MTAB_6701_obj, 
                               metadata = E_MTAB_6701_metadata)
# check the metadata
E_MTAB_6701_obj[[]]
unique(E_MTAB_6701_obj$final_cluster)
unique(E_MTAB_6701_obj$location)

E_MTAB_6701_placenta <- subset(E_MTAB_6701_obj, subset= location %in% "Placenta" )
E_MTAB_6701_placenta[[]]
unique(E_MTAB_6701_placenta$final_cluster)
unique(E_MTAB_6701_placenta$annotation)
unique(E_MTAB_6701_placenta$location)

dput(unique(E_MTAB_6701_placenta$annotation))
dput(unique(E_MTAB_6701_placenta$final_cluster))

saveRDS(E_MTAB_6701_obj, "E_MTAB_6701_obj.rds")
E_MTAB_6678_obj <- readRDS("E_MTAB_6678_obj.rds")
E_MTAB_6701_obj <- readRDS("E_MTAB_6701_obj.rds")

E_MTAB_6678_obj[[]]
E_MTAB_6701_obj[[]]
rownames(E_MTAB_6701_obj)
unique(E_MTAB_6701_obj$location)

placenta_sub <- subset(E_MTAB_6701_obj, subset = location %in% "Placenta")
saveRDS(placenta_sub, "placenta.only.rds")

placenta_merge <- merge(x= E_MTAB_6678_obj,
                        y= E_MTAB_6701_obj,
                        add.cell.ids = c("E_MTAB_6678","E_MTAB_6701"))


placenta_merge <- NormalizeData(placenta_merge)
placenta_merge <- FindVariableFeatures(placenta_merge)
placenta_merge <- ScaleData(placenta_merge)
placenta_merge <- RunPCA(placenta_merge)
ElbowPlot(placenta_merge)
placenta_merge <- RunUMAP(placenta_merge, dims=1:20)

DimPlot(placenta_merge)

placenta_int <- IntegrateLayers(object = placenta_merge, 
                            method = CCAIntegration, 
                            orig.reduction = "pca", 
                            new.reduction = "integrated.cca",
                            verbose = FALSE)

placenta_int[["RNA"]] <- JoinLayers(placenta_int[["RNA"]])

placenta_int <- RunUMAP(placenta_int, 
                        reduction = "integrated.cca", 
                        dims = 1:20)
saveRDS(placenta_int, "vento_placenta_int.rds")
placenta_int[[]]

colors <- c("#00cc99","#6f00ff","#ff9933","#17becf","#cc00cc","#b19cd9","#f7b6d2","#0070ff",
            "#ff2052","#e78ac3","#00bfff","#be0032","#b5a642","#9acd32","#ff00ff","#008000",
            "#ff55a3","#e49b0f","#ff6347","#915c83","#738678","#c49c94","#1f77b4","#ff7f0e", 
            "#c54b8c","#ffbf00","#9467bd","#eee600","#e377c2","#7f7f7f","#bcbd22","#aec7e8",
            "#fdee00","#a7fc00","#c5b0d5","#e3dac9","#c7c7c7","#dbdb8d","#9edae5","#ad494a",
            "#8c6d31","#bd9e39","#e7ba52","#7b4173","#843c39","#2ca02c","#d62728","#1B9E77",
            "#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666","#fddde6",
            "#86608e","#49796b","#fefe22","#a52a2a","#e41a1c","#377eb8","#4daf4a","#984ea3",
            "#ff7f00","#ffff33","#a65628","#f781bf","#999999","#66c2a5","#fc8d62","#8da0cb", 
            "#f0dc82","#a6d854","#ffd92f","#e5c494","#b3b3b3","#8dd3c7","#c9a0dc","#0073cf",
            "#ff9896","#f6a600","#98777b","#8c564b","#273be2","#cc0000","#0070ff","#fb8072")

Idents(placenta_int) <- "annotation"
DimPlot(placenta_int, cols = colors, alpha = 0.4, label = T, repel = T)
placenta_int <- FindNeighbors(placenta_int)
rownames(placenta_int)
FeaturePlot(placenta_int,
            features = "GBP2-ENSG00000162645")

features <- c('ERVFRD-1-ENSG00000244476',
              'ERVMER34-1-ENSG00000226887',
              'ERVW-1-ENSG00000242950',
              'FCGR2A-ENSG00000143226',
              'FCGR3A-ENSG00000203747',
              'FCGRT-ENSG00000104870',
              'FURIN-ENSG00000140564',
              'GBP2-ENSG00000162645',
              'GBP5-ENSG00000154451',
              'IFITM3-ENSG00000142089',
              'MARCH8-ENSG00000165406',
              'PCSK1-ENSG00000175426',
              'PCSK2-ENSG00000125851',
              'PCSK5-ENSG00000099139',
              'PCSK6-ENSG00000140479',
              'PCSK7-ENSG00000160613',
              'SERINC5-ENSG00000164300')

FeaturePlot_scCustom(placenta_int,
                     features = features)


# Define the output directory
output_dir <- "feature_plots"
dir.create(output_dir, showWarnings = FALSE)  # Create directory if not exist

# Loop through each feature and save the plot
for (feature in features) {
  p <- FeaturePlot_scCustom(placenta_int, features = feature)
  
  ggsave(filename = file.path(output_dir, paste0(feature, ".png")),
         plot = p,
         width = 6, height = 6, dpi = 300)  # High resolution (300 dpi)
}

saveRDS(placenta_int, "vento_placenta_int.rds")

placenta_int <- readRDS("vento_placenta_int.rds")

placenta_sub <- subset(placenta_int, subset= location == "Placenta")
Idents(placenta_sub) <- "annotation" 



genes <- c(
           'PCSK1-ENSG00000175426',
           'PCSK2-ENSG00000125851',
           'FURIN-ENSG00000140564',
           'PCSK5-ENSG00000099139',
           'PCSK6-ENSG00000140479',
           'PCSK7-ENSG00000160613',
           
           'GBP2-ENSG00000162645',
           'GBP5-ENSG00000154451',
           'IFITM3-ENSG00000142089',
           'MARCH8-ENSG00000165406',
           'SERINC5-ENSG00000164300',
           'MARCH8-ENSG00000165406') 

average_expression <- AverageExpression(placenta_sub, 
                                        features = genes, 
                                        assays = "RNA", 
                                        slot = "data")


expression_data <- average_expression$RNA

# Reorder rows to match the order in `genes`
expression_data <- expression_data[genes, , drop = FALSE]

head(expression_data)
expression_data[1:4,1:4]

# Scale expression data (column-wise)
scaled_expression_data <- t(scale(t(expression_data)))

library(viridis)
library(colorRamps)
library(circlize)
library(RColorBrewer)

# Define viridis color mapping
col_fun <- colorRamp2(c(min(scaled_expression_data), 0, max(scaled_expression_data)), 
                      viridis(3))

col_fun <- colorRamp2(c(min(scaled_expression_data), 0, max(scaled_expression_data)), 
                      rev(brewer.pal(3, "Spectral")))


# Generate the heatmap with viridis colors
heatmap <- Heatmap(as.matrix(scaled_expression_data), 
                   name = "Scaled Expression",
                   col = col_fun,  # Apply viridis colors
                   show_row_dend = TRUE,
                   cluster_columns = FALSE,
                   cluster_rows = FALSE,
                   column_names_side  = "top",
                   row_names_gp = gpar(fontsize = 10))

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# Your gene list, split into two groups
genes_group1 <- c('PCSK1-ENSG00000175426', 'PCSK2-ENSG00000125851', 'FURIN-ENSG00000140564',
                  'PCSK5-ENSG00000099139', 'PCSK6-ENSG00000140479', 'PCSK7-ENSG00000160613')

genes_group2 <- c('GBP2-ENSG00000162645', 'GBP5-ENSG00000154451', 'IFITM3-ENSG00000142089',
                  'MARCH8-ENSG00000165406', 'SERINC5-ENSG00000164300', 'MARCH8-ENSG00000165406')

# Define the groups for splitting
gene_groups <- c(rep("group1", length(genes_group1)), 
                 rep("group2", length(genes_group2)))

# Extract and order expression data accordingly
expression_data <- expression_data[c(genes_group1, genes_group2), , drop = FALSE]

# Scale expression data (column-wise)
scaled_expression_data <- t(scale(t(expression_data)))

col_fun <- colorRamp2(c(min(scaled_expression_data), 0, max(scaled_expression_data)), 
                      viridis(3))

col_fun <- colorRamp2(c(-4, 0, 4), viridis(3)) 

# Generate the heatmap with row splits
Heatmap(as.matrix(scaled_expression_data), 
        name = "Expression",
        col = col_fun,
        cluster_rows = FALSE,  # Keep the order as defined
        cluster_columns = FALSE,  
        row_split = (gene_groups),  # Use split for row separation
        column_names_side  = "top",
        show_column_dend = FALSE,
        row_names_gp = gpar(fontsize = 10))

# complete
