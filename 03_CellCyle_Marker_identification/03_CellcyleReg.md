
### Import Libraries

```{R}
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

```

Firstly we will normalise the data with respect to readcount to account for sequencing depth variation across the cells. `sctransform` function of seurat help to achieve both normalisation and Variance stabilisation. Before we regress out the variation due to sequencing it may be a good idea to check if there is celle cycle based variation in out dataset.

In order to find the Cell-cycle based variation we will need a rough normalisation of data and this we will achieve by dividing by total cellcounts per cell and then taking the natural log.  

```{R}
# Normalize the counts
wt1Filtered<- NormalizeData(wt1Filtered)
```

First we will fetch the Cell cycle specific gene IDs.  Code and explanation text below is taken from **https://hbctraining.github.io/scRNA-seq/lessons/cell_cycle_scoring.html** .  The suggested link also have cellcycle gens for Zebrafish, Drosophilla and Humans at link **https://github.com/hbc/tinyatlas/tree/master/cell_cycle**. 

*All of the cell cycle genes are Ensembl IDs, but our gene IDs are the gene names. To score the genes in our count matrix for cell cycle, we need to obtain the gene names for the cell cycle genes.We can use annotation databases to acquire these IDs. While there are many different options, including BioMart, AnnotationDBI, and AnnotationHub. We will use the AnnotationHub R package to query Ensembl using the ensembldb R package.*

```{R}
# Download cell cycle genes for organism at https://github.com/hbc/tinyatlas/tree/master/cell_cycle. Read it in with:

cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv") 
cell_cycle_genes <- read.csv(text = cc_file)


#BiocManager::install("AnnotationHub")
#BiocManager::install("ensembldb")

library(AnnotationHub)
# Connect to AnnotationHub
am <- AnnotationHub()

# Access the Ensembl database for organism
amDb <- query(am, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- amDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

# Download the appropriate Ensembldb database
edb <- am[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)
```

*Now we can use these annotations to get the corresponding gene names for the Ensembl IDs of the cell cycle genes.*

```{R}

cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "S") %>%
        pull("gene_name")
        
# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "G2/M") %>%
        pull("gene_name")
```

*Taking the gene names for the cell cycle genes we can score each cell based which stage of the cell cycle it is most likely to be in. By default, the PCA is run only using the most variable features identified previously. The output of the PCA returns the correlated gene sets associated with the different principal components (PCs).*

```{R}
wt1Filtered <- CellCycleScoring(wt1Filtered,
                                   g2m.features = g2m_genes,
                                   s.features = s_genes)


# View cell cycle scores and phases assigned to cells                   
View(wt1Filtered@meta.data)  

# Identify the most variable genes
wt1Filtered <- FindVariableFeatures(wt1Filtered, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
wt1Filtered <- ScaleData(wt1Filtered)		     

# Perform PCA
wt1Filtered <- RunPCA(wt1Filtered)

# Plot the PCA colored by cell cycle phase

DimPlot(wt1Filtered,
        reduction = "pca",
        group.by= "Phase")



DimPlot(wt1Filtered,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
```


Here we do not see large differences due to cell cycle phase. Based on this plot, we would not regress out the variation due to cell cycle.

#### SCTransform

Now we will use the `sctranform` function will carry out a proper normalisation accounting for sequencing depth of each cell (or nUMIs).  Since we donot have to account for Cell cycle we will regress out the variation due to mitichondrial genes using the `vars.to regress` argument of the function.

Dueing this process the R object may grow in size and we have to increase the size of allowable R object size.  This parameter changes with the size of input dataset. default setting is **500 X 1024 ^ 2 = 500 Mb**. Beloww we are raising that to ~2GB.
```{R}
options(future.globals.maxSize = 8000 * 1024^2)
```

```{R}
wt1Filtered<-SCTransform(wt1Filtered, vars.to.regress = c("percent.mt"))
```
The next step will to dimentionality reduction where we would be first calculating Principle components and then will use the relevant PC for graphical clustering using UMAP.  Here we will use the information from all the key PCA and then plot in two dimensions using UMAP (Uniform Manifold Approximation and Projection).
Lets start with calculating PCAs and then we can plot the PCAs

```{R}
# Run PCA
seurat_integrated <- RunPCA(object = seurat_integrated)

# Plot the PCA, we can do PCA a plot of multiple samples by using `split.by = "sample"` argument
PCAPlot(seurat_integrated)

```


The next step is to run UMAP using defined number of PCA's serving as dimensions in the function.  However, we need to decide how many PCAs we need to use for UPA.  There are couple of different ways.  One is to do `ElbowPlot` of PCAs

```{R}
ElbowPlot(wt1Filtered,ndim=50)
```

In the graph we see an elbow in the curve around PC12-PC16 that accounts for the majority of the variance.  After that the components do carry the variance of dataset but their contribution can be low. Generally where the graph flattens out is used a refrence for choosing PC for UMAP.

There are two other ways to choose PCs.  One is visual where you plot a heatmap of different PC's. This will help us visualize which PCs are relevant and explain significant variation in the dataset.
The code below will plot first 10 PCAs using top 500 cells.

```{R}
DimHeatmap(wt1Filtered, dims = 1:9, cells = 500, balanced = TRUE)
```
Other way to decide is to do `ElbowPlot` with cummulative Standard deviation vs percent Std deviation.  Choose a PC where 
1. The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
2.The point where the percent change in variation between the consecutive PCs is less than 0.1%.

```{R}
# Determine percent of variation associated with each PC
pct <- wt1Filtered[["pca"]]@stdev / sum(seurat_integrated[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

 #Create a dataframe with values
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))

# Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text(size=2) + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()


```

Now lets run UMAP with first 40PC's in default mode before we try playing with parameters and do the whole process step by step.

```{R}
wt1Filtered<-RunUMAP(wt1Filtered, 
                dims = 1:40,
		        reduction = "pca")

DimPlot(wt1Filtered)   

```

We can also look for gene expression across these cluster using `FeaturePlot` function.

```{R}
FeaturePlot(wt1Filtered,features=c("Nrxn3", "Slco1a4", "Npy","Gpc5" ), cols = c("lightgray","red"))

```

#### Cluster Cells

This are the steps that goes in UMAP clustering if we do step by step.  Seurat uses Graph Based Clustering  which embeds cells in a Graph Structure using K-nearest neighbor Graph (KNN) to cluster cell with similar gene expression.

This Occur in multiple steps.
(1) Step1: In the first step we will calculate the nearest neighbours. 
```{R}
wt1Filtered <- FindNeighbors(object = wt1Filtered, 
                                dims = 1:41)
```

(2) Step2: Find clusters.  Here we will use `FindClusters` function to find the cluster.  An important parameter of the function is `resolution`. The resolution parameter sets the "granularity" of the downstream clustering and will need to be optimized for every individual experiment. For datasets of 3,000 - 5,000 cells, the resolution set between 0.4-1.4 generally yields good clustering. Increased resolution values lead to a greater number of clusters, which is often required for larger datasets.

```{R}
# Determine the K-nearest neighbor graph
wt1Filtered  <- FindNeighbors(object = wt1Filtered , 
                                dims = 1:41)
                                
# Determine the clusters for various resolutions                                
wt1Filtered <- FindClusters(object = wt1Filtered ,
                               resolution = c(0.25,0.4, 0.6, 0.8, 1.0, 1.4))

```

The assignment of each cell to cluster can be identified from the metadata of seurat object `wt1Filtered`.

```{R}
wt1Filtered@meta.data%>% View()
```

We have tried couple of different resolutions. we can have tried couple of different resolutions.  lets start with 0.6. We will assign identities to cluster based on the the resolution. By changing `integrated_snn_res.NN` the identities of the cluster can be changed.
```{R}
# Assign identity of clusters
Idents(object = wt1Filtered) <- "SCT_snn_res.0.25"

# Plot the UMAP
DimPlot(wt1Filtered,
        reduction = "umap",
        label = TRUE,
        label.size = 6)


```
Here we can use the `FeaturePlot` to explore the status of additional metrics from Metadata.

```{R}
metrics<-c("nCount_RNA","nFeature_RNA","percent.mt","G2M.Score", "G2M.Score")

FeaturePlot(wt1Filtered, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
```

The feature plot can also help to look at the expression of specific genes. This can help us identify cluster of specific cell type.

```{R}
FeaturePlot(wt1Filtered, 
            reduction = "umap", 
            features = c("Nrxn3","Npy"),
            pt.size = 0.4, 
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE)
```

#### Marker Identification
The marker identification of cluster can be done at 3 levels
1.  Identify all markers for each cluster.  This compares any given cluster with all the other clusters and outputs the genes that are differentialy expressed in the cluster. *Useful for identifying unknown clusters and improving confidence in hypothesized cell types.* 

2. Identification of conserved markers of the cluster: This is applicable when we have multiple samples. This analysis identify genes that are differentially expressed/present within each sample first, and then reports those genes that are conserved in the cluster across all conditions. These genes can help to figure out the identity for the cluster.*Useful with more than one condition to identify cell type markers that are conserved across conditions.*

3. Marker identification between specific clusters: this analysis explores differentially expressed genes between specific clusters.(*Useful for determining differences in gene expression between clusters that appear to be representing the same celltype (i.e with markers that are similar) from the above analyses.*)

4. Marker identification between specific cell types. In this scenario the 2 group of cells may not belong to distinct clusters and in that case we may have to identify the barcode of cells and then use them to find markers.

The function to find markers is `FindAllMarkers` and have 3 important arguments

`logfc.threshold` :  Minimum log2FC for average expression incluster/average expression in other clusters combined [default=0.25]. 
**WARNING** 
(a)If the gene is expressed in only a smal fraction of cell and does not meeet threshold that will be missed.
(b) return a lot of metabolic/ribosomal genes due to slight differences in metabolic output by different cell types, which are not as useful to distinguish cell type identities

`min.diff.pct` : minimum percent difference between the percent of cells expressing the gene in the cluster and the percent of cells expressing gene in all other clusters combined.

`min.pct`: only test genes that are detected in a minimum fraction of cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1.

These parameters will setup the stringency of identifying the markers. The `only.pos` is set to `TRUE` to report only positive changes.

```{R}

# Find markers for every cluster compared to all remaining cells, report only the positive ones
markers <- FindAllMarkers(object = wt1Filtered, 
                          only.pos = TRUE,
                          logfc.threshold = 0.5)                     

```

Run `View(markers)` to check identified markers.  If we have multiple samples we can execute `FindConservedMarkers()` to identify markers that are conserved between clusters when compared to other cluster in the sample. For more details **https://github.com/hbctraining/scRNA-seq/blob/master/lessons/09_merged_SC_marker_identification.md**

For a good marker we want hight pct.1 and low pct.2 (check the `markers` object )


The code below can help to identify the Cell barcodes if cells that belown to a givenn cluster.
```{R}
colnames(wt1Filtered)[wt1Filtered$seurat_clusters==0]
```


Lets find the **top10** markers of each cluster.
```{R}
markers %>% group_by(cluster) %>% top_n(10,avg_log2FC)

```
Now we can use the `FeaturePlot` function to check the expression of markers in the cluster.    

```{R}
top_genes<-markers %>% group_by(cluster) %>% top_n(2,avg_log2FC)%>%pull(gene)
FeaturePlot(wt1Filtered,
            features = top_genes[1:3],
            col=c("lightgray","red"),
            pt.size=0.1)
```

**Other Visualisations**
```{R}
VlnPlot(wt1Filtered,features = top_genes[1:3], pt.size=0.5)

RidgePlot(wt1Filtered, features = top_genes[1:3])
```

Finding differentially expressed genes between 2 clusters.

```{R}
dge1vs2<-FindMarkers(wt1Filtered,ident.1 = 1,ident.2 = 2, only.pos = FALSE, logfc.threshold = 0)

# Ratio will be cluster1/cluster2
# only.pos = FALSE : To get all genes for GSEA
# logfc.threshold = 0 : All the fold changes
# Multiple cluster can be specified by providing vector of cluster
# e.g. ident.2 = c(2,3,6,7)
```

We can Rename the identities of the cell, once we have identified cell types and then can use `DimPlot` to visualise clusters with new label.

```{R}
# Rename all identities (Names are arbitrary can be replaced by real cell identities.)
wt1Filtered <- RenameIdents(object = wt1Filtered, 
                               "0" = "Rango",
                               "1" = "Mexico",
                               "2" = "DieHard",
                               "3" = "SpiderMan",
                               "4" = "DarkNight",
                               "5" = "IronMan",
                               "6" = "Gaurdians",
                               "7" = "Galaxy",
                               "8" = "Deadpool",
                               "9" = "Wolverine",
                               "10" = "Mask",
                               "11" = "Eraser",
                               "12" = "Recall",
                               "13" = "StarTrek",
                               "14" = "Starwars",
                               "15" = "Jedai",
			                   "16" = "OK",
			                   "17" = "Summer", 
			                   "18" = "Spring", 
			                   "19" = "GotIT!")

# Plot the UMAP
DimPlot(object = wt1Filtered, 
        reduction = "umap", 
        label = TRUE,
        label.size = 3,
        repel = TRUE)

```
The `DotPlot` function can show the expression of markers along different clusters.  Lets pick up top 2 markers from each cluster and plot them in `DotPlot`.

```{R}
#Getting Top 2 features
top2m<-markers %>% group_by(cluster) %>% top_n(2,avg_log2FC)%>%pull(gene)
top2m
# Dotplot

DotPlot(wt1Filtered,features=rev(top2m), cols=c("blue","red"), dot.scale=2)+RotatedAxis()

```

### Custom cell type analysis.
In this scenari we are interested in certain cell types which donot form distinct clusters but are scattered over multiple clusters.  In this case we will look into 3 cell types 

- Celltype 1 : Cells expressing Dcx, Foxp2, Meis2 (not Sp8, Pax6)
- Celltype 2 : Cells expressing Dcx, Sp8 (not Foxp2, Meis2, Pax6)
- Celltype 3 : Cells expressing Dcx, Pax6 (not Foxp2, Meis2,Sp8)

Before we look for markers specific for each cell type lets investigate the expression of each of these markers in cells.
```{R}
FeaturePlot(wt1Filtered, features = c("Foxp2","Sp8","Pax6","Dcx","Meis2"))

VlnPlot(wt1Filtered, features = c("Foxp2","Sp8","Pax6","Dcx"))
```
In order to account for expression values we have to assign RNA measurements as default assay for Seurat object.

```{R}
DefaultAssay(wt1Filtered) <- "RNA"
```
Now lets investigate the expression of each marker.
```{R}
p1<-FeaturePlot(wt1Filtered,features=c("Dcx" ),min.cutoff=0, cols = c("lightgray","red"))

p2<-FeaturePlot(wt1Filtered,features=c("Foxp2","Meis2" ),min.cutoff=0, cols = c("lightgray","red"))

p3<-FeaturePlot(wt1Filtered,features=c("Sp8" ),min.cutoff=0, cols = c("lightgray","red"))

p4<-FeaturePlot(wt1Filtered,features=c("Pax6" ),min.cutoff=0, cols = c("lightgray","red"))

VlnPlot(wt1Filtered, features = c("Foxp2","Sp8","Pax6","Dcx","Meis2"))

plot_grid(p1, p2, p3, p4)
```
We observed that there may be overlap of expression of markers. To identify the cells expression our criteria we have to do a bit of clever coding shown below. First all 3 celltype express Dcx, so lets isolate all the cells expressing Dcx.

```{R}
dcx<-subset(wt1Filtered,subset= (Dcx >0))

FeaturePlot(dcx,features=c("Dcx" ))

DefaultAssay(dcx) <- "RNA"
```

Now lets find the cell types that are expressing genes of our interest.

```{R}
# Find cell barcodes expressing each marker.
foxp<-rownames((subset(dcx,subset= (rna_Foxp2 >0)))@meta.data)
sp<-rownames((subset(dcx,subset= (rna_Sp8 >0)))@meta.data)
pax<-rownames((subset(dcx,subset= (rna_Pax6 >0)))@meta.data)
meis<-rownames((subset(dcx,subset= (rna_Meis2 >0)))@meta.data)

# List of all barcodes in dcx object
dcx_barcodes<-rownames(dcx@meta.data)

# BIT OF WIERD CODING
# The idea here is to generate combinatorial expression of markers
f1<-ifelse(dcx_barcodes %in% foxp,1,0)
s2<-ifelse(dcx_barcodes %in% sp,2,0)
p4<-ifelse(dcx_barcodes %in% pax,4,0)
m8<-ifelse(dcx_barcodes %in% meis,8,0)
marker_expression<-f1+s2+p4+m8

dcx[["marker_expression"]]<-marker_expression

head(dcx@meta.data)
tail(dcx@meta.data)

DimPlot(dcx,label=TRUE)

# Cell barcodes expressing Dcx, Foxp2 and Meis2
foxp2_9<-rownames(dcx@meta.data[dcx@meta.data$marker_expression2==1,])
length(foxp2_9)

# Cell barcodes expressing Dcx and Sp8
sp8_2<-rownames(dcx@meta.data[dcx@meta.data$marker_expression2==2,])
length(sp8_2)

# Cell barcodes expressing Dcx and Foxp2
pax6_4<-rownames(dcx@meta.data[dcx@meta.data$marker_expression2==4,])
length(pax6_4)

```
Now that we have barcode of cells expressing the markers we are interested in lets calculate the Markers(DE genes).  The code is below to calculate DE genes and saving them as CSV files.

```{R}
wtcombined_foxp2_9_vs_sp8_2.markers <- FindMarkers(wt.combined, ident.1=foxp2_9, ident.2 = sp8_2, min.pct = 0.25)

wtcombined_foxp2_9_vs_pax6_4.markers <- FindMarkers(wt.combined, ident.1=foxp2_9, ident.2 = pax6_4, min.pct = 0.25)

wtcombined_sp8_2_vs_pax6_4.markers <- FindMarkers(wt.combined, ident.1=sp8_2, ident.2 = pax6_4, min.pct = 0.25)
```







