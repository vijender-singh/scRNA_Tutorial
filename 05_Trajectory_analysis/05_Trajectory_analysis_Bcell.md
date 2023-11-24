
## Trajectory Analysis (pseudotime analysis)

**Acknowledgement:** The script has been reproduced from Khusbu Patel's github repo https://github.com/kpatel427 .  I have added some additional codes to explain some parts a bit better. 

The dataset we will be using for demonstration purpose is published in the article (link below) studying transcriptomic landscape of human blood cells using singlecell technology.  The study has 7551 human blood cells covering 32 cell types from 21 healthy donors but for the sake of demonstration today we will not be using the entire data set we are just choosing a group of cells a group of B cells and its progenitors which are around 1400 cells so we will be building a trajectory using these group of cells. So the goal of today's analysis is to construct a trajectory to order cells in pseudotime and finally find the genes that change expression as cells progress along a directory. 

[Data](http://scrna.sklehabc.com/): http://scrna.sklehabc.com/

[Alternate Data Link](https://drive.google.com/drive/folders/1eJZIaWOMgKllFquIt9-W60SiuHLDny26?usp=sharing)

[Publication](https://academic.oup.com/nsr/article/8/3/nwaa180/5896476?login=false): https://academic.oup.com/nsr/article/8/3/nwaa180/5896476?login=false

[Monocle3 tutorial](https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/): https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/

R package collection for [Trajectory Inference](https://dynverse.org/): https://dynverse.org/

[Publication]( https://www.biorxiv.org/content/10.1101/276907v1.full.pdf) comparing various Trajectory Inference methods: https://www.biorxiv.org/content/10.1101/276907v1.full.pdf


so for the demonstration today we will be using following packages
```{R}
set.seed(1234)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)
```
First download the [Data](http://scrna.sklehabc.com/) from the Data or [Alternate Data Link](https://drive.google.com/drive/folders/1eJZIaWOMgKllFquIt9-W60SiuHLDny26?usp=sharing).  

```{R}
# read in data from Data link
markers <- read.delim('ABC_Marker.txt', header = T) # gene metadata
metadata <- read.delim('ABC_Meta.txt', header = T) # cell metadata
expr <- read.delim('ABC_umi_matrix_7551_cells.csv', header = T, sep = ',') # expression matrix
head(markers)
head(metadata)
head(expr[1:10.1:10])

```
Data downloaded from alternate source can be importde in R using the code below. This will import exactly the 3 objects `markers`,`metadata` and `expr` that are created in the above code.
```{R}
# For data from alternate link can be loaded in R studio as
load("Monocle3_Tutorials_files.rdata") 
head(markers)
head(metadata)
head(expr[1:10.1:10])
```
In our analysis we will be doing intial filtering and QC in **seurat** and then we will move to using Monocle for Trajectory analysis.

**Did you notice something different about the expression matrix `expr`?**

After having a closer look, we see that it is transposed expression matrix where Cell ids are Rows and genes as Coloumns.  However, for using this in seurat we have to transpose back so that Cells ids are columns and genes are rows. Thats what we do in first line of code `t(expr)`. 

```{R}
# create seurat object
expr.t <- t(expr)
dim(expr.t)
head(expr.t[1:10.1:10])
```

After that we are creating a seurat object as we have been doing in other tutorials. Following that lets visualise the metdata of this seurat object.
```{R}
seu.obj <- CreateSeuratObject(counts = expr.t)
View(seu.obj@meta.data)
```
For the metadata we can see that it has the number of counts and number of features for each of the cell ids. We also need to add additional metadata here which includes the annotation and a bunch of information that we have in the cell metadata so let us add the cell metadata information to this metadata. To achieve that we can merge these sarat objects metadata with our cell metadata by merging the two by column. For seurat metadata we will use the `row.names` and `cell_id` for cell metadata table. And then we will assign the newly merged metadata table to seurats metadata (over writing the old one). While doing that we will assign `Row.names` column as rownames using `column_to_rownames` function.
```{R}
seu.obj@meta.data <- merge(seu.obj@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id')
View(seu.obj@meta.data)

seu.obj@meta.data <- seu.obj@meta.data %>% 
  column_to_rownames(var = 'Row.names')
```
Next, lets calculate Mitochondrial gene fraction in each cell.  Following that we filter the seurat object with our usual parameters of minimum features, minimum UMIs and max mitochondrial content.
```{R}
seu.obj$mitopercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')

seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 &
                    nFeature_RNA > 500 &
                    mitopercent < 10)

# Lets check how many cells we begin with and how many of them we filtered out.

seu.obj

seu.obj.filtered
```
The next step is to subset only B-cells from all the cell types and name this new seurat object as `b.seu`. Before that lets see what all annotation/cell types are there in dataset.

```{R}
# subset my seurat object - B cells

unique(seu.obj.filtered@meta.data$population)

table(seu.obj.filtered@meta.data$population)
```
The celltype in the population are 
- sp : Stem cell progenitor
- t :  T- cells
- mo: Monocytes
- nk: Natural Killer cells
- e : Erythrocytes
- b : B-cells
- n : Neutrophills

Lets subser seurat object for B cells.  But lets check what is the default identity of the cells is set to.  
```{R}
Idents(seu.obj.filtered)
```
Ah !  it is set to some other identity, we want to change it to the population identity so that we can filter these cells based on the identity `b`. This can be achived by `Idents(seurat.object) <- NewIdentity` function.
```{R}
Idents(seu.obj.filtered) <- seu.obj.filtered$population
# Lets check if the identities are now set to population.
Idents(seu.obj.filtered)
```
Now we can filter B cells as,

```{R}
b.seu <- subset(seu.obj.filtered, idents = "b")
b.seu
```
The filtered B-cells have additional types or identities or class assigned to them. Lets check that
```{R}
unique(b.seu@meta.data$redefined_cluster)
table(b.seu@meta.data$redefined_cluster)
```
We can see that in this B-cell population we have cells ranging in various stages of differentiation/development from progenior B-ceel to terminal Plasma or Memory B-cell.  This is ideal dataset to perform Trajectory analysis.

We will now do the regular steps we do for dimentionality reduction `NormalizeData`, `FindVariableFeatures`, `ScaleData`, `RunPCA`, `FindNeighbors`, `FindClusters` and finally `RunUMAP`.

**NOTE:** *While finding clusters using `FindClusters` we are using resolution 0.9. I want to point out that trajectory analysis is heavily dependent on the topology of the data and the topology is determined by the clustering. Hence, it is very important that you cluster your data using the optimal resolution making sure that each cell type groups together as distinct clusters and dis-similar groups of cells do not get clubbed together in one cluster.  This might be important as there might be transition states of cells which could be easily grouped together if we're using a lower resolution so make sure you're using the optimal resolution to cluster your data*

```{R}
# pre-processing using seurat
b.seu <- NormalizeData(b.seu)
b.seu <- FindVariableFeatures(b.seu)
b.seu <- ScaleData(b.seu)
b.seu <- RunPCA(b.seu)
b.seu <- FindNeighbors(b.seu, dims = 1:30)
b.seu <- FindClusters(b.seu, resolution = 0.9)
b.seu <- RunUMAP(b.seu, dims = 1:30, n.neighbors = 50)
```

```{R}
a1 <- DimPlot(b.seu, reduction = 'umap', group.by = 'redefined_cluster', label = T)
a2 <- DimPlot(b.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2
```
**Note:** *Here the clustering resolution has done a fair job in terms of separating or grouping the major cell types into distinct clusters. Although we do see some amount of overlaps between some cell types but I think some major cell types like probe cycling preview memory plasma cells are grouped into separate clusters. So the recommendation is that you play around with resolutions when you're clustering your data to ensure that you do not over cluster or under cluster your data.  Be sure that the major cell types in your data do group together into distinct clusters.*

**How to find the optimal resolution?** 

*It's an iterative process so we run our clustering steps with various resolutions and then visualize our data as a UMAP. We can also take a look at the feature plot looking at the expression of certain markers that we know are expressed in certain cell types to make sure that these cells are separated out in distinct clusters.  We iterate the process with different resolutions until we do see distinct clusters for distinct cell types.*

Here I think my clustering resolution has done a fair job although it's not ideal because we do see some clusters that are not technically distinct clusters. 

**NOTE:** If we are doing your own analysis, we have to spend some time annotating different cluster into cell types or stages as per our need.

At this point we have everythoing to move to Monocle3 for performing Trajectory analysis.  Monocle requires an object of `cell_data_set` class.  we need to convert the Seurat object into an object of `cell_data_set` class. Seurat wrappers has a function called `as.cell_data_set` to achieve that. `as.cell_data_set` object is based on single cell experiment class the way we retrieve certain information is different than how you do it in Seurat.  We will extract some of these information below

```{R}
# Convert to cell_data_set object 

cds <- as.cell_data_set(b.seu)
# Lets visualise tyhe object
cds

# to get cell metadata we use colData
colData(cds)

# to gene metdata we use fdata 
fData(cds)
# Gene metadat is a dataframe with 19813 rows with rownames but ) column.
# We can extract the first 10 gene names
rownames(fData(cds))[1:10]

# Here we will create a column gene_short_name that holds the gene names which are present as row names.
# Why? : Is it easy to extract these names if they are column id dataframe then merely rownames.

fData(cds)$gene_short_name <- rownames(fData(cds))
fData(cds)

# to get counts as sparse matrix
counts(cds)
```
**- RECAP**

We have converted the Seurat object to an object of `cell_data_set` class, we also learned how to retrieve information from `cell_data_set` object.  

**-Whats Next?**

Now to build a trajectory we want monocle3 to use the clustering information previously performed in Seurat. For that we need to retrieve clustering information and store it inside appropriate locations in the `cell_data_set` object in monocle. Also monocle has functions to cluster cells where it not only determines the clusters but also partitions and these partitions are nothing but these are superclusters
So we need to change three things, 

**- 1** Need to add the partition information i.e, Assign Paritions (Super Cluster)*. 

**- 2**  Add the clustering information from Seurat.

**- 3**  Include UMAP cell embeddings information which are UMAP coordinates.

**Assign Paritions (Super Cluster)**

`cds` object has multiple slots which are different than what Seurat has. In cds the clustering information is in `cds@clusters$UMAP${cluster_result,partition,clusters}`. For partition we want to assign all cells in same partition/super cluster as we are expecting only one trajectory. If we are expecting multiple trahectories or branches we will have more superclusters. For this we have to create a supercluster labelled 1 (all cells belong to this cluster). This will be done by creating a named factor vector whose length is same as number of cells and value is 1.

```{R}
# Cluster cells (using clustering info from seurat's UMAP)
# let's use the clustering information have

# Creating Named Vector
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

# Assigning named factor vector to partition slot
cds@clusters$UMAP$partitions <- reacreate.partition
cds@clusters$UMAP$partitions 
```
**Add the clustering information from Seurat.**

Next we extract the cluster information from seurat and add to the slot in `cds` object at `cds@clusters$UMAP$clusters`

```{R}
# Assign the cluster info
#Cluster information 
b.seu@active.ident

# Extracting cluster information and saving to list_cluster variable.
list_cluster <- b.seu@active.ident

# Assignint cluster information to correct slot in cds.
cds@clusters$UMAP$clusters <- list_cluster
```
**Adding UMAP embedding information**

In this step we will extract the UMPA embedding information, i.e, their coordinate location in 2D UMAP plot. The slot location for this is `cds@int_colData@listData$reducedDims${UMAP,PCA}` where emebedding information for both PCA and UMAP are stored.  We assign the UMPA information.

```{R}
# Assign UMAP coordinate - cell embeddings
# UMAP embeddings in Seurat
b.seu@reductions$umap@cell.embeddings

# Assigning that in cds slot
cds@int_colData@listData$reducedDims$UMAP <- b.seu@reductions$umap@cell.embeddings
```
Similar to `DimPlot` of Seurat Monocle3 has `plot_cells` function to visualise the cells.  Here we are visualising the cells twice once by labelling with `cluster` (generated in seurat) and another one with `redefined_cluster` (came with the dataset).   

```{R}
# plot

cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE, 
           group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
           color_cells_by = "redefined_cluster",
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names
```

Now the next step in the analysis is to learn the trajectory graph using the `learn_graph` and we will set the `use_partition = FALSE` as we donot expect multiple trajectories or branches. If this is set to true then it will use partitions and learn a disjoint trajectory for each partition but since we know that all our cells um follow one trajectory we will set this to false.

```{R}
#  Learn trajectory graph 
cds <- learn_graph(cds, use_partition = FALSE)
```
Lets visualise the Data
```{R}
plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)

```
Looking at the plot here we can see that a trajectory has been learnt which is going from one end of the cells to the other and also all of the cells are connected by one trajectory.  we do not see any disjointed trajectories in addition to that we also see some branch points here.  A possible explanation could be that when we were grouping/clustering cells we could not distinctly separate out all the cell types into separate clusters and hence we could be seeing these branch points.  Also there could be cells from the later stage or the earlier stage being present in cluster and hence the branching. This would make more sense once we know which cells are at the beginning of the trajectory which are at the later stage by ordering these cells in pseudotime so the next step would be to order the cells in pseudotime. 

Pseudotime is estimated using `order_cells` function on `cds` and specifying the reduction method (UMAP here) and also specifying the  `root_cells`.  `root_cells` are the precursor cell list in fporm of a vector from where the rest of the lineages developed/orginated/differentoiated. This has to be the prior knowledge.  Here we know that pro-B cells are the start and these cells are located in cluster 5.      

```{R}
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 5]))
```

Lets visualise our cells based on pseudotime and color the cells based on pseudotime they are associated with.
```{R}
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)

```
In plot we can see that the pro-B cells which we provided to the monocle as root cells have the lowest pseudo time and as cell progresses along the trajectory they have higher and higher pseudo time with plasam cells having the highest pseudotimes.  Now let us visualize the pseudotime value for each type of cell. We can get pseudotime values using the `pseudotime` function and it gives a pseudotime value for each cell. We can store this information in variable called `monocle3_pseudo` time of `cds` object. This isnformation will be present as a part of the `colData` attribute.

```{R}
# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
colData(cds)
```
To explore the pseudotime values in each of the `redefined_cluster`, first we can create a dataframe of all the metadata associated with `colData(cds)` and then visualise it as box plot.

```{R}
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, aes(monocle3_pseudotime, fill = redefined_cluster)) +
  geom_boxplot()
```
The box plot shows range of pseudotime associated with each cell type.  Lets order the boxplots a bit better.
```
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(redefined_cluster, monocle3_pseudotime, median), fill = redefined_cluster)) +
  geom_boxplot()
```
The cells are arranged with median value of pseudotime. We can see that the pro-B cells have the lowest pseudotime and the plasma cells have the highest median pseudotime.  This gives us a sense of trajectory as to what are the various stages B-cell go through starting from pro-B cells to plasma cells. This trajectory or this lineage also corroborates with the literature that was experimentally validated.  So this trajectory analysis is now helpful to look at genes that change their expression starting from pro-B cells and going through different stages to ultimately becoming the plasma B cells.  It would be interesting to identify genes that are differentially regulated or activated or repressed during the earlier or the later stages of this trajectory.  

To identify genes that change expression as cells progress along a trajectory we will basically be using a function called `graph_test` that tests the genes for differential expression based on the low dimension embeddings.  It uses a test statistics that tells whether the cells at nearby positions on a trajectory will have similar or dissimilar expression levels for the gene that is being tested.  So let's run the graph test function with the first parameter is `cds` object for neighbor_graph we provide principal_graph and we provide four cores. 

```{R}
# Finding genes that change as a function of pseudotime.  Code is same as suggested in the manual 
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)
view(deg_bcells)
```

Lets identify the genes that have highest statistical significance of expression changes during the pseudotime tramnsition between stages.

```{R}
deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()
```
Now that we have list of significant genes that changes their expression.  Lets visualise them using `Seurat`'s `FeaturePlot` function.
```{R}
FeaturePlot(b.seu, features = c('E2F2', 'STMN1', 'CD52'))
```
**visualizing pseudotime in seurat**
We can add the pseudotime information to Seurat and visualise our new information on Seurat object.
```{R}
b.seu$pseudotime <- pseudotime(cds)
FeaturePlot(b.seu, features = "pseudotime")
```
The Seurat FeaturePlot is now showing the psudotime attribute of the cells.  Lets add cell cluster identities and replot the `FeaturePlot`.
```{R}
Idents(b.seu) <- b.seu$redefined_cluster
FeaturePlot(b.seu, features = "pseudotime", label = T)
```
