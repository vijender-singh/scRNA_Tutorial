## scRNA_Tutorial

**NOTE :** This single-cell tutorial is dervied from different sources and I would like to thank the individuals maintaining the resources and making them available for wider research community. The sources for each tutorial are acknowledged in each individual tutorial section.

### Tutorial Layout
The tutorial covers Data download from SRA and cellranger preprocessing to generate count matrix.  The count matrix then can be imported in Seurat to create a Seurat object. Besides the data download and Cellranger tutorial covers three aspects of single cell analysis
- Generating count Matrix: [Data download and Cellranger](https://github.com/vijender-singh/scRNA_Tutorial/blob/main/01_DataDownload/01_DataDownload_processing.md) (fastq to count matrix)
- Single Sample scRNA analysis.
    - [ScRNA Quality check](https://github.com/vijender-singh/scRNA_Tutorial/blob/main/02_DataQC_Filtering/02_SingleCellTutorial_filtered.md)
    - [scRNA Cellcycle stage evaluation and analysis](https://github.com/vijender-singh/scRNA_Tutorial/blob/main/03_CellCyle_Marker_identification/03_CellcyleReg.md)
- [Case-Control scRNA sample analysis](https://github.com/vijender-singh/scRNA_Tutorial/blob/main/04_Case_control_integrateData/04_DataIntegration_CaseControl.md)
- [Trajectory analysis](https://github.com/vijender-singh/scRNA_Tutorial/blob/main/05_Trajectory_analysis/05_Trajectory_analysis_Bcell.md)

The markdown file for each section of the tutorial is in each folder.
