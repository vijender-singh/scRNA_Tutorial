# CountMatrix QC and filteration.

## Source Declaration  
The tutorial is mostly derived from the example and tutorial material from Seurat website. The dataset used here is taken from Tepe et al., [2018] study where they used scRNA to understand cellular heterogeneity and functionality during neurogenesis in olfactory bulb. The accession number for the raw single-cell RNA sequencing data reported in this paper is GEO: GSE121891.  We will be using the WT1 sample for understanding basic processing and analysis of scRNA dataset.

Lets load the required packages.
```{R}
#BiocManager::install("seurat")
suppressMessages(require("Seurat", quietly = TRUE))
suppressMessages(require(dplyr, quietly = TRUE))
suppressMessages(require(patchwork, quietly = TRUE))
suppressMessages(require(data.table, quietly = TRUE))
suppressMessages(require(tidyverse, quietly = TRUE))
```

## Data import
The output of `cellranger count` function will produce  `barcodes.tsv.gz`, `features.tsv.gz` and `matrix.mtx.gz` files for each 10X sample.  These files can be imported into R using `Seurat`'s `Read10x` function.  However this function doesn't allow any prefix to be added to these file names to differentiate between output of multiple samples. One way to come around is to have separate folders for each sample with respective three `cellranger count` output files (`barcodes.tsv.gz`, `features.tsv.gz` and `matrix.mtx.gz`).  I personally like to have a prefix before the files, so I modified the orginal function to `Read10xprefix` that allows to specify the prefix of the files.  The `Read10xprefix`  is present in `scRNA_custom.R` script and to use the function we will `source` this file.

Before importing these files in R lets explore what each file contains. 
```{R}
head(fread("./WT1_filtered_feature_bc_matrix/features.tsv.gz"))
```
`fread` is a function from  `data.table` package that can read `.gz` compressed files and then head prints the first 6 rows.  So here we see that the first column contains mouse Ensemble geneID, second column have their common names and lastely third column indicates nature of the feature.

```{R}
head(fread("./WT1_raw_feature_bc_matrix/WT1_barcodes.tsv.gz"),10)
```
This contains single column that have cellular barcodeID's and each barcode represent single cell. Ignore the `-1`

```{R}
#head(fread("./WT1_raw_feature_bc_matrix/WT1_matrix.mtx.gz")[1:10.1:10])

```
This file contains the counts for each transcript in each cell.  Since the data from single cell is zero inflated  as large number of genes have no counts for them hence the value is zero.  Storing the zeros take too much space so we store the data in sparse matrix format. We will take all these thre files and convert this sparse matrix into dense matrix format that will appear as a nice matrix of 2 dimension. 

This is shown in the image below.  
![image](./denseMatrix_ScRNA2.png)

## Creating `Seurat` object
First we will create a dense matrix and then will use it to create Seurat object.

```{R}
# This will allow to import prefixed files.
source("scRNA_custom.R")
wt1.data<-Read10xprefix(data.dir="./WT1_raw_feature_bc_matrix", prefix="WT1_")
class(wt1.data)
dim(wt1.data)
head(wt1.data[1:5,1:5])
```
So we have created a dense matrix (also refered as **countmatrix**)  of 32285 rows (# of genes) and 536398 columns (cells).

```{R}
#Lets look at CreateSeuratObject function
?CreateSeuratObject

#Lets apply the basic filtering
# min.cells =3   # Only considergenes which are observed in minimum 3cells

wt1<-CreateSeuratObject(counts = wt1.data, project = "WT1", min.cells = 3, min.features = 200)
wt1
```


# Quality Control and Filtering

**Doublets** :  Doublet are the situation when a single droplet is catching more than one cell.  The proportion of these cells will vary from the tissue/cell types.  The problem that doublets pose is to have a unique transcription profile  which sometimes indicates of cells that may be transitioning between two states/types. This could be misleading and can lead to wild goose chase in analysis.

**How are doublets identified?**

They general metric for identification of a doublet is detection of large number of genes than an average. This seems to be a reasonable argument but that it will falsely identify the transitioning cells as doublets.There are tools out there that can help with identifying doublets like `Scrublet`.  We will not be using in this tutorial.

Before moving into actual analysis we would like to refine out countmatrix by filtering out low quality cells. Low quality cells refer to cell being stressed, broken, dead or droplets with multiple cells. Since they misrepresent the information on a cell we will exclude them from analysis.  A detailed account on cell filtering and QC is published by Ilicic et. al.[2016]. 
Quoting from the article 

"*There is an extensive literature on the relationship between mtDNA, mitochondrially localized proteins, and cell death . Upregulation of RNA levels of mtDNA/RNA in broken cells suggests losses in cytoplasmic content. In such situation where cell membrane is broken, cytoplasmic RNA will be lost, but RNAs enclosed in the mitochondria will be retained*"

Such broken cell can be identified with higher mitochondrial RNA content.
We will apply filter to clean following issues
1. Broken, empty droplets are characterized by detection of very few genes.
2. Droplets with multiple cells will show higher gene content. This will also show increase in detection of unique genes.
3. Dying or ruptured cells show higher mitochondrial contamination. The mitochondrial QC metrics can calculated with `PercentageFeatureSet()` function. This function calculated counts originating from a subset of features. We will specify these features as mitochondrial genes by  specifying `MT-` prefix. In some gene annotation this is in lowercase as `mt-`, so change accordingly.

In nut shell we are using 3 (1 optional) criteria to filter cells
- How many uniques genes detected in the cell?
- How many uniques molecules detected in each cell?
- What fraction of counts are assigned to Mitochondrial genes?
- How many genes are detected per UMI? This will indicate the complexity of the data (optional).

Lets calculate fraction of mitochondrial reads.

```{R}

wt1[["percent.mt"]] <- PercentageFeatureSet(wt1, pattern = "^mt-")
VlnPlot(wt1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0.25,ncol = 3)

## Lets have a look at metadata

head(wt1@meta.data)
```

Metadata has following columns

`orig.ident`: Sample identity if known, otherwise assigned the `project` variable.

`nCount_RNA`: Number of UMIs per cell

`nFeature_RNA`: Number of genes detected per cell

We can make a better sense of this data if we plot Mitochondrial fraction `percent.mt` against `nFeature_RNA` and `nCount_RNA`. This can be achieved either by `qplot` function of `tidyverse` or  `FeatureScatter` function of Seurat.


```{R}
plot1 <- FeatureScatter(wt1,feature1="nCount_RNA",feature2="percent.mt")

plot2 <- FeatureScatter(wt1,feature1="nCount_RNA",feature2="nFeature_RNA")

plot1 + plot2
 
```
Theres is one cell that has outlier counts and making the plot squeeze towards the y axis.  Lets remove it and replot the data.

```{R}
### Identify the cellID based on max nCount_RNA
cellID<-rownames(wt1@meta.data[which((wt1@meta.data[["nCount_RNA"]])==max(wt1@meta.data[["nCount_RNA"]])),])

### Remove cell from the dataset and plot with ggplot2
library(ggplot2)
wt1Filtered <- wt1[,!colnames(wt1) %in% cellID]
plot1 <- ggplot(wt1Filtered@meta.data,aes(x=log10(nCount_RNA),y=percent.mt))+
         geom_point(col="red",alpha=0.2)+
         xlab("UMI/cell")+ylab("percent Mitochondrial")
plot2 <- ggplot(wt1Filtered@meta.data,aes(x=log10(nCount_RNA),y=nFeature_RNA))+
         geom_point(col="red",alpha=0.2)+
         xlab("UMI/cell")+ylab("Genes/cell")
plot3<- ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA,y=percent.mt))+
         geom_point(col="red",alpha=0.2)+
         xlab("Genes/cell")+ylab("percent Mitochondrial")+
         geom_vline(xintercept = 250)

plot1+plot2+plot3
```
### **Assessing the Quality metrics**

**UMI Counts per cell**  

This is a good metric to use for filtering generally we can extract good information if cells are above 500 UMI's.  UMI's between 500-1000 is preferred as it can provide usable information. However a higher read depth sequencing is preferred.

```{R}
ggplot(wt1Filtered@meta.data,aes(x=nCount_RNA,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 500)
```

**Genes Detected per Cell**

The next parameter to look at is number of genes detected per cell.  This cutoff is generally around 200-300 range but can be made less or more stringent but never less than 200. Cells with less the 250 genes are possibly lysed/broken and this can be seen in the above graphs (Genes/cell vs percent Mitochondrial). The graph below is far from ideal, a high quality data will have a single peak that encapsulates most of the cells in the assay.  We see a second peak around 1500 this may indicate a distinct population of cells which are larger in size and may have higher gene content.

```{R}
ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA)) + 
  	geom_density(alpha = 0.2,col="#5BB300",fill="#5BB300") + 
  	scale_x_log10() + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 250)
```

**UMI vs Genes Detected**

Another important metric to look at is UMI vs genes detected.  We expect to see an increase in UMIs/Cell with increase in Genes/Cell.  We will also include the mitochondrial content of each cell. As observed that the cells with lower gene content tend to have higher mitochondrial fraction.  Our UMI and Gene filter will be unable to remove these cells so we have to filter the cells based on the mitochondrial content of the cell.


```{R}
ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA,y=nCount_RNA,col=percent.mt))+
  geom_point()+
  scale_x_log10() + 
  scale_y_log10() +
  stat_smooth(method=lm,col="yellow")+
  geom_vline(xintercept=300)+geom_hline(yintercept=500)+geom_vline(xintercept=250,col="red")
```

**Mitochondrial Fraction as density plot**

As described before higher MT gene content is sign of dying, lysed or poor quality cells. The Cutoff can be stringent as 5% but max should be no higher that 20% unless we are expecting a higher % of mitochondrial content in cells (cardiac muscles ~5000 mitochondria/cell, 100,000-300,000 mitochondria per ovum)


```{R}
ggplot(wt1Filtered@meta.data,aes(x=percent.mt,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 20, col="red")+	geom_vline(xintercept = 12.5, col="blue")

```

**Complexity**

This is measured as Genes per UMI. A sample will be highly complex if every single UMI represent a distinct gene, in that case the score/ratio is 1.  However, in a typical cell this ratio is expected around 0.8. A cell with low complexity indicates a very smaller set of genes are expressed e.g like RBCs.

```{R}
wt1Filtered[["GenesPerUMI"]]<-(log10(wt1Filtered@meta.data$nFeature_RNA)/log10(wt1Filtered@meta.data$nCount_RNA))

ggplot(wt1Filtered@meta.data,aes(x=GenesPerUMI,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	ylab("density") +
  	geom_vline(xintercept = 0.8, col="blue")+
    xlab("log10(Genes per UMI)")

```
**Reads per cell**

This is another important metric to look at when comparing 2 samples.  Ideally both samples should peak at roughly same location between 10,000-100,000reads per cell.  The information for this plot requires information from mapping steps to be stored and used.


**Filtering**

We have explore the different aspects of the single cell data.  Based on the biology of the cell type the filtering parameters has to be set up.  Here we will apply the following filter

nUMI > 500

nGene > 300

log10GenesperUMI >0.8

Mitochondrial percentage < 20

```{R}
wt1Filtered<-subset(wt1Filtered,
            subset=(nFeature_RNA >=200) & (nCount_RNA >=500) & (percent.mt < 20) & (GenesPerUMI > 0.80))

wt1Filtered2<-subset(wt1,
            subset=(nFeature_RNA >=200) & (nCount_RNA >=500) & (percent.mt < 20))# & (GenesPerUMI > 0.80))

```

**Gene Level Filtering**

In our filtered cells there will be genes that have zero counts. So we will apply 2 level filter
1) Remove genes that has zero counts across all cells
2) Remove genes that are expressed in less than 10 cells ( Careful, if we expect very few cells to express a particular gene than this filter has to be adjusted accordingly.)

```{R}
# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = wt1Filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
wt1Filtered <- CreateSeuratObject(filtered_counts, project="WT1",meta.data = wt1Filtered@meta.data)

```


### RE-ASSESS QC Metrics

```{R}

wt1Filtered[["percent.mt"]] <- PercentageFeatureSet(wt1, pattern = "^mt-")
VlnPlot(wt1Filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0.25,ncol = 3)

## Lets ahve a look at metadata
```

### UMI, Genes and Mitochondrial content

```{R}
plot1 <- ggplot(wt1Filtered@meta.data,aes(x=log10(nCount_RNA),y=percent.mt))+
         geom_point(col="red",alpha=0.2)+
         xlab("UMI/cell")+ylab("percent Mitochondrial")
plot2 <- ggplot(wt1Filtered@meta.data,aes(x=log10(nCount_RNA),y=nFeature_RNA))+
         geom_point(col="red",alpha=0.2)+
         xlab("UMI/cell")+ylab("Genes/cell")
plot3<- ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA,y=percent.mt))+
         geom_point(col="red",alpha=0.2)+
         xlab("Genes/cell")+ylab("percent Mitochondrial")+
         geom_vline(xintercept = 250)

plot1+plot2+plot3
```

### Distribution of Genes Detected

```{R}
ggplot(wt1Filtered@meta.data,aes(x=nCount_RNA,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 500)
```
```{R}
ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA)) + 
  	geom_density(alpha = 0.2,col="#5BB300",fill="#5BB300") + 
  	scale_x_log10() + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 250)
```

### Genes Detected vs UMI Content per cell

```{R}
ggplot(wt1Filtered@meta.data,aes(x=nFeature_RNA,y=nCount_RNA,col=percent.mt))+
  geom_point()+
  scale_x_log10() + 
  scale_y_log10() +
  stat_smooth(method=lm,col="yellow")+
  geom_vline(xintercept=300)+geom_hline(yintercept=500)+geom_vline(xintercept=250,col="red")
```
### Percent Mitochondrial Content

```{R}
ggplot(wt1Filtered@meta.data,aes(x=percent.mt,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	ylab("log 10 Cell density") +
  	geom_vline(xintercept = 20, col="red")+	geom_vline(xintercept = 12.5, col="blue")

```

### Genes per UMI

```{R}
wt1Filtered[["GenesPerUMI"]]<-(log10(wt1Filtered@meta.data$nFeature_RNA)/log10(wt1Filtered@meta.data$nCount_RNA))

ggplot(wt1Filtered@meta.data,aes(x=GenesPerUMI,fill="red")) + 
  	geom_density(alpha = 0.2) + 
  	ylab("density") +
  	geom_vline(xintercept = 0.8, col="blue")+
    xlab("log10(Genes per UMI)")

```

Save the filtered object for downstream analysis
```{R}
save(wt1Filtered, file="WT1_filtered.rdata")
```



Publications:

1.Tepe B, Hill MC, Pekarek BT, Hunt PJ, Martin TJ, Martin JF, Arenkiel BR. Single-Cell RNA-Seq of Mouse Olfactory Bulb Reveals Cellular Heterogeneity and Activity-Dependent Molecular Census of Adult-Born Neurons. Cell Rep. 2018 Dec 4;25(10):2689-2703

2.Ilicic, Tomislav et al. “Classification of low quality cells from single-cell RNA-seq data.” Genome biology vol. 17 29. 17 Feb. 2016, doi:10.1186/s13059-016-0888-1

### Acknowledgement:
A large chunk of the code is either taken from Seurat website or from  https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html with few or sometimes no modification.  I personally would like to thank the owners of both the resources for provinding such a great material on single cell.    

