# create the second object

# read in data

ast_data <- Read10X("/Users/rhalenathomas/Documents/Data/scRNAseq/AST23_BrainComm/CellRangerOuts/AST23/raw_feature_bc_matrix")
dim(ast_data)


#Look at a small part of the data
ast_data[1:5, 1:5]

#Look at the distribution of the number of UMIs per cell
colSums(ast_data) %>% summary
#Look at the distribution of the number of genes per cell  
colSums(ast_data > 0) %>% summary

#Remove barcodes with less than 100 genes detected (you can select a different value here)
# without 
ast_data <- ast_data[, colSums(ast_data > 0)> 100]
dim(ast_data)

# create Seurat object
ast_seu <- CreateSeuratObject(ast_data, project = "SNCA", min.cells = 3)
# look at the object dimensions
ast_seu

# look at the distribution of total counts of RNA across cells
ast_seu$nCount_RNA %>% summary

# look at the distribution of unique RNA transcripts across cells

ast_seu$nFeature_RNA %>% summary

#Calculate the percentage of RNA encoded mitochondrial genes from the mitochondrial DNA
ast_seu <- PercentageFeatureSet(ast_seu, pattern = "^MT-", col.name = "percent.MT")
ast_seu$percent.MT %>% summary

VlnPlot(ast_seu, features = "percent.MT", pt.size = 0.001)

# filter object with the same settings as with the control
ast_seu_ft <- subset(ast_seu, nCount_RNA < 60000 & nFeature_RNA > 500 & 
                       percent.MT < 20)


# see the results
VlnPlot(ast_seu, features = c("percent.MT", "nCount_RNA", "nFeature_RNA"), pt.size = 0.001)

# check how many cells we have
dim(ast_seu_ft)

# Normalize data (log normalization) and select genes with variable expression across cells --------------------------------------

ast_seu_ft <- NormalizeData(ast_seu_ft, normalization.method = "LogNormalize", scale.factor = 10000)

#Check out the effect of normalization
GetAssayData(ast_seu_ft, assay = "RNA", slot = "data") %>% expm1 %>% colSums %>% head
GetAssayData(ast_seu_ft, assay = "RNA", slot = "counts") %>% colSums %>% head

# PCA analysis
# Find Variable features
# our selection method is vst
ast_seu_ft <- FindVariableFeatures(ast_seu_ft, selection.method = "vst", nfeatures = 2000)

var  <- VariableFeatures(ast_seu_ft)
VariableFeaturePlot(ast_seu_ft)
var[1:10]

#Scaling is recommended before PCA, as otherwise highly expressed genes will have a disproportionate effect on the PC composition

# we are also regressing MT genes to remove them from the PCA
ast_seu_ft <- ScaleData(ast_seu_ft, vars.to.regress = "percent.MT")
ast_seu_ft@assays$RNA@scale.data %>% dim

#Linear dimensionality reduction
#Choosing the number of PCs can depend on how many cells you have
ast_seu_ft <- RunPCA(ast_seu_ft, assay = "RNA", npcs = 30)

PCAPlot(ast_seu_ft)

#Assess how many PCs capture most of the information in the data 
ElbowPlot(ast_seu_ft, ndims = 30)

# visualize with UMAP
ast_seu_ft <- RunUMAP(ast_seu_ft, dims = 1:18)
UMAPPlot(ast_seu_ft)

# Remove doublets with Doublet finder
sweep.res.ast <- paramSweep_v3(ast_seu_ft, PCs = 1:18, sct = FALSE)

#We do not have the "ground truth" regarding doublets, such from from genotype data for pooled samples 
#We summarize the performance of the range of pN=pK parameters we tested
sweep.stats_ast <- summarizeSweep(sweep.res.ast, GT = FALSE)

#Here the "best" pK for the data is chosen based on a metric determined by the DoubletFinder developers
#Which performs best in datasets where the ground truth is known
bcmvn_ast <- find.pK(sweep.stats_ast)
ggplot(bcmvn_ast, aes(x = pK, y = BCmetric, group = "Sweep")) + geom_point() + geom_line() + 
  theme(axis.text.x = element_text(angle = 90))


#We will pick pK = 0.15
#We are not going to use our clustering information to estimate "homotypic" doublets
#We are simply going to use an expected doublet formation rate of 2.5% based on the number of starting cells loaded

nExp_poi <- round(0.025*nrow(ast_seu_ft@meta.data))
ast_seu_ft <- doubletFinder_v3(ast_seu_ft, PCs = 1:18, 
                               pN = 0.25, pK = 0.15, 
                               nExp = nExp_poi, 
                               reuse.pANN = FALSE, sct = FALSE)

#Here we update the Seurat object version just in case the one returned by DoubletFinder is an older version
ast_seu_ft <- UpdateSeuratObject(ast_seu_ft)

#Visualize and assess the cells called as probable doublets
UMAPPlot(ast_seu_ft, group.by = "DF.classifications_0.25_0.15_30")

# table of doublets and signlets
ast_seu_ft$DF.classifications_0.25_0.15_30 %>% table

# visualize the features in doublets and singlets
VlnPlot(ast_seu_ft, features = c("nCount_RNA", "nFeature_RNA", 
                                 "percent.MT", "pANN_0.25_0.15_30"), 
        group.by = "DF.classifications_0.25_0.15_30", pt.size = 0.001, ncol = 2, log = TRUE)

# we have identified doublets but we haven't removed them at this point
# remove the doublets
# we can do this by subset
# first we need to set the active meta data slot to the doublet identification
Idents(ast_seu_ft) <- "DF.classifications_0.25_0.15_30"
# we select only the singlet cells
ast_seu_ft2 <- subset(ast_seu_ft, idents = "Singlet")

