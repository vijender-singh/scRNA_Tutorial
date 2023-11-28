# Running Cell Ranger to generate CountMatrix

A count matrix is a matrix whose rows represent each genes and columns represent individual cell, and cells of the matrix will hold the expression of a gene in a given cell. We will use the `cellranger count` function of [**`cellranger`**](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/tutorials/gex-analysis-nature-publication#GEXIntro). In order to run `cellranger count` function we need to have fastq files either from the sequencer or downloaded form SRA.

## Scenario 1: Use of Published datasets
 We will look into how to download files from SRA as sometimes to do a comparative analysis of your data with published data you mayneed to download fastq files.

### STEPS DOWNLOADING FASTQ FILES

#### Step1 : Create a file with accession numbers

Create a file on cluster `accessionlist.txt` (can be named differently) that contain the accession number of SRA runs using nano.  For tutorial purpose we will download the data published in [Tepe et al., 2018, Cell Reports 25, 2689–2703](https://www.sciencedirect.com/science/article/pii/S2211124718317972?via%3Dihub).  Here the authors characterized the cellular heterogeneity and neuron development in the olfactory bulb of mouse.  We will use only the WT samples.  We will start by creating the file `accessionlist.txt` listing SRA run files corresponding to **WT1** (SRR8126390-SRR8126393) and **WT2** (SRR8126394-SRR8126397)samples.

```
SRR8126390
SRR8126391
SRR8126392
SRR8126393
SRR8126394
SRR8126395
SRR8126396
SRR8126397
```
#### Step2 : Download the data with a SRAtoolkit
Next we will be dowloading all the SRA files listed in the `accessionlist.txt` with `fastq-dump` function of SRAtoolkit.  We will submit an array job to download the files on the cluster.

```
#!/bin/bash
#SBATCH --job-name=datadownload
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 6
#SBATCH --mem=18G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --array=[0-7]
#SBATCH -o ./logs/%x_%A_%a.out

hostname

# Load SRAtoolkit/2.11.3 module
module load sratoolkit/2.11.3

# Create an array holding all the SRA accession numbers list
SRAlist=($(cat accession.txt))

# SLURM_ARRAY_TASK_ID will help us pick one specific accession number from the array.
# That specific accession number is assigned to runID variable.
runID=${SRAlist[$SLURM_ARRAY_TASK_ID]}

# Next three lines of code create a TMPDIR so that if the default TMPDIR has issues it
#doesnot affect our run.
pathTMP=`pwd`
mkdir -p ${pathTMP}/tmp
export TMPDIR=${pathTMP}/tmp

# fastq-dump to download the files.
fastq-dump --split-files --gzip ${runID}

# --split-files  : Create seperate files for R1 (labelled _1) and R2 (labelled _2) reads
# --gzip         : Compress the files
```

Once the downloads are complete we will have multiple runs, 8 to be precise, for each sample (**WT1** and **WT2**).  It is better we merge these 8 files to create 2 files for each sample, one for R1 (fwd) and another for R2 (Rev) reads. The code below is to achieve the merging of files. The merged files should follow `Sample_S1_L00X_R1_001.fastq.gz` naming convention for their use in downstream steps. **L00X** indicates sequencer lane information, since we donot have this information on this data we are assigning all reads to `L001`.

```
#!/bin/bash
#SBATCH -J merge
#SBATCH -p general
#SBATCH -q general
#SBATCH --mem=50G
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -o ./logs/%x_%j.out

cat SRR8126394_1* SRR8126395_1* SRR8126396_1* SRR8126397_1* > WT2_S1_L001_R1_001.fastq.gz
cat SRR8126394_2* SRR8126395_2* SRR8126396_2* SRR8126397_2* > WT2_S1_L001_R2_001.fastq.gz

cat SRR8126390_1* SRR8126391_1* SRR8126392_1* SRR8126393_1* > WT1_S1_L001_R1_001.fastq.gz
cat SRR8126390_2* SRR8126391_2* SRR8126392_2* SRR8126393_2* > WT1_S1_L001_R2_001.fastq.gz
```
At this stage we have fastq files that we can directly use in `cellranger count` function to create count matrix.

## Scenario 2: Recieved data from Sequencer
Sequencing facility may deliever our data in the form of fastq files or raw data. In later case, from raw data we have to generate the fastq files.  This is done by using `cellranger mkfastq` function. Following the files you will see when you download the data.

```
Files   Files.metadata  Properties  Samples
```
Of the four directories the `Files` directory is the one we are interested in and it contains the intensities that have to be converted to fastq files. IN order to do that we have to create a sample csv file (cellranger-tiny-bcl-simple-1.2.0.csv). The layout of the file is shown below.
```
Lane,Sample,Index
1,ARID1A,SI-GA-A11
1,WT2,SI-GA-B11
2,ARID1A,SI-GA-C11
2,WT2,SI-GA-D11
```
Here there are 2 samples (2nd column WT2 and ARID1A) from a project that were sequenced together on 2 lanes, lane 1 and 2 (first column) using the Indexes mentioned in column3. The information lannes and Indexes has to be obtaioned from the sequencing facility.

That is all we need to generate fastq files.  Now lets looks at the script with `cellranger mkfastq` code.
```
#!/bin/bash
#SBATCH -J CrMkFq
#SBATCH -c 8
#SBATCH -p general
#SBATCH --qos=general
#SBATCH --mem=100G
#SBATCH -o ./logs/%x_%j.out

module load cellranger/2.1.1
module load bcl2fastq/2.20

# Here we are giving project an arbitrary name "project1"

cellranger mkfastq --id=project1 \
--run=/path/to/Files \
--csv=cellranger-tiny-bcl-simple-1.2.0.csv \
--output-dir=/path/to/project1/fastq \
--localcores=8 --localmem=100
```
Once we exceute this script, 2 output directories will be created with the fastq files at the location specified in `--output-dir` option.  So the path will look something like
```
/path/to/project1/fastq/XNXXN/ARID1A
/path/to/project1/fastq/XNXXN/WT2
```
The XNXXN is a directory name with (cahracters and numbers in uppercase). In my understanding this is completely random. The above path will have all the fastq files.

Now that we have the fastq files, lets create the countsmatrix. In order to map reads from our fastq files we need reference transcriptome and that has to be indexed before we can use it.  10X genomics maintain indexed reference transcriptome for Mouse and Human on their [website](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).  The website has instructions of downloading them. Lets assume that we have downloaded the mouse indexes and they are located at `/path/to/resources/refdata-gex-mm10-2020-A/`

Below is the script that can be used to create the count matrix.

```
#!/bin/bash
#SBATCH -p xeon
#SBATCH -q general
#SBATCH -c 8
#SBATCH --mem=80G
#SBATCH -o ./logs/counts-%j.out

cd ../02_count

module load cellranger/7.0.0

cellranger count --id=WT1 \
                   --transcriptome=/path/to/resources/refdata-gex-mm10-2020-A/ \
                   --fastqs=/path/to/merged/fastq \
                   --sample=WT1 \
                   --localcores=8 \
                   --localmem=80

```
 The output will create a WT1 outpu tdirectory and this directory will have an `outs` directory, that will have following contents
 ```
analysis                       metrics_summary.csv       possorted_genome_bam.bam.bai
filtered_feature_bc_matrix     molecule_info.h5          raw_feature_bc_matrix
filtered_feature_bc_matrix.h5  possorted_genome_bam.bam  raw_feature_bc_matrix.h5
 ```

 `analysis` : This directory will have stats on PCA, t-SNE, UMPA plots in default mode.

 `filtered_feature_bc_matrix` : This directory contains the count matrix with filtered data.

 `raw_feature_bc_matrix` : This directory contains the count matrix with raw data.

Following will be the files in `filtered_feature_bc_matrix` or `raw_feature_bc_matrix`.

```
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```
`barcodes.tsv.gz` : Information on Cell barcodes.
`features.tsv.gz `: Information on features (genes/transcripts) that are quantified.
`matrix.mtx.gz`   :  This contains actual countdata in sparse matrix format.

So we now have our starting material, count matrix to be further analysed.  The downstream analysis will be using `bioconductor` and `CRAN` packages in `R`.  So the best otion is to download this data and then analyse in R studio.  The use of cluster resources end here.  However if one wants to continue on cluster the rest of analysis can be done on cluster as well.
