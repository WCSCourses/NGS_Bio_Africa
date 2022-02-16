# RNA-Seq Expression Analysis: _Bacteria_
## 1.1 Introduction
RNA sequencing (RNA-Seq) is a high-throughput method used to profile the transcriptome, quantify gene expression and discover novel RNA molecules. This tutorial uses RNA sequencing of from two _**Mycobacterium tuberculosis**_ isolates (each with 3 replicates - 6 total samples) to walk you through transcriptome pseudo-alignment and quantification and how to profile transcriptomic differences between experimental conditions (**Sensitive** vs. **Resistant**) by identifying differentially expressed genes.

For an introduction to RNA-Seq principles and best practices see:  
> **A survey of best practices for RNA-Seq data analysis.** Ana Conesa, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, Daniel J. Gaffney, Laura L. Elo, Xuegong Zhang and Ali Mortazavi Genome Biol. 2016 Jan 26;17:13 doi:10.1186/s13059-016-0881-8

## 1.2. Learning Outcomes
By the end of this practical, you should be familiar with:  
<ul>
    <li> Using <code>salmon</code> to quantify transcript abundance 
    <li> Knowing how to create a study design file
    <li> Loading the required files into <code>R</code>
    <li> Creating a <code>DESeq2</code> data set
    <li> Running <code>DESeq</code>
    <li> Visualising the relationship between samples 
    <li> Creating a <code>results</code> object
    <li> Ranking genes by significance
    <li> Visualisation of differential expression results
</ul>

Focus on:  
<ul>
    <li> Looking at what is in the files you are given, and the files you create 
    <li> How the information from different files relates and comes together 
    <li> Testing different inputs to commands to see what happens
</ul>

## 1.3. Practical Outline
<p align="center">
		<img src="https://github.com/WCSCourses/NGS_Bio_Africa/blob/main/images/H3ABioNet_Logo%20(1).png" style="width:100%">
		<b>Workflow Overview.</b> A schematic representation of each step in this practical exercise.
</p>

The main steps in all differential expression analyses are:  
<ul>
    <li> <b>Quantification</b> - How many reads come from each gene? (<code>Salmon</code>) 
    <li> <b>Normalisation</b> - Dealing with biases in the data (<code>DESeq2</code>)
    <li> <b>Differential expression analysis</b> (<code>DESeq2</code>)
    <li> <b>Visualisation and reporting</b> (<code>R libraries</code>)
</ul>

>**Additional Resources:**  
><a href="https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html">RNA-seq workflow: gene-level exploratory analysis and differential expression</a>

## 1.4. Authors
This tutorial was developed by Jon Ambler, Phelelani Mpangase and Nyasha Chambwe, based in part, from materials from Victoria Offord and Adam Reid.

## 1.5. Prerequisites
This tutorial assumes that you have the following software or packages and their dependencies installed on your computer. The software or packages used in this tutorial may be updated from time to time so, we have also given you the version which was used when writing the tutorial. Where necessary, instructions for how to download additional analysis tools are given in the relevant section.

<div align="center">
    <table>
        <tr>
            <th width="20%"> Package </th>
			<th width="60%"> Link for download/installation instructions </th>
			<th width="20%"> Version tested</th>
		</tr>
		<tr>
			<td> Trimmmomatic </td>
			<td> <a href="https://github.com/usadellab/Trimmomatic">https://github.com/usadellab/Trimmomatic</a> </td>
			<td> 0.39 </td>
		</tr>
		<tr>
			<td> Salmon </td>
			<td> <a href="https://salmon.readthedocs.io/en/latest/index.html">https://salmon.readthedocs.io/en/latest/index.html</a> </td>
			<td> 1.4.0</td>
		</tr>
		<tr>
			<td> R </td>
			<td> <a href="https://www.r-project.org/">https://www.r-project.org/</a> </td>
			<td> 4.0.2 </td>
		</tr>
		<tr>
			<td> tximport </td>
			<td> <a href="https://bioconductor.org/packages/release/bioc/html/tximport.html">https://bioconductor.org/packages/release/bioc/html/tximport.html</a> </td>
			<td> 1.20.0 </td>
		</tr>
		<tr>
			<td> GenomicFeatures </td>
			<td> <a href="https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html">https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html</a> </td>
			<td> 1.44.0 </td>
		</tr>
		<tr>
			<td> DESeq2</td>
			<td> <a href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html">https://bioconductor.org/packages/release/bioc/html/DESeq2.html</a> </td>
			<td> 1.32 </td>
		</tr>
		<tr>
			<td> pheatmap </td>
			<td> <a href="https://cran.r-project.org/web/packages/pheatmap/index.html">https://cran.r-project.org/web/packages/pheatmap/index.html</a> </td>
			<td> 1.0.12 </td>
		</tr>
	</table>
</div>

## 1.6. Setup
### 1.6.1. Create Practical Directory
Navigate to the module folder on the Virtual Machine (VM) using the following command:  
```bash
cd /home/manager/course_data/rna_seq_pathogen
```

Create a working directory for the practical:  
```bash
mkdir practical
```

Create a copy of the practial dataset to maintain the data integrity of the original dataset:  
```bash
 cp /home/manager/course_data/rna_seq_pathogen/data/bacterial/* /home/manager/course_data/rna_seq_pathogen/practical
```

Copy and paste the command above as a single line at your command line window.  
>Note: if you make a mistake at any point during this tutorial, you can reset by deleting the practical folder and restarting from this section.

Move to the practical directory:  
```bash
cd practical
```

### 1.6.2. Download Supplemental Datasets
**Internet Access Required**

#### 1.6.2.1. Reference Transcriptome File
Download the GFF formatted version of the annotation file. Unfortunately due to the lack of standardization in bioinformatics, we need both the <a href="https://mblab.wustl.edu/GTF22.html">>GTF</a> file and the <a href="https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md">GFF</a> file for analysis (See the differences and similarities between these two file format specifications <a href="https://www.ensembl.org/info/website/upload/gff.html?redirect=no">here</a>). Salmon uses the GTF file and the GenomicFeatures library in R needs a GFF file. Both files contain the same annotation information, they are just formatted differently. No tools are able to reliably convert one to the other because even GFF files are not formatted consistently.  
```bash
wget https://www.dropbox.com/s/4yjgbmy3dyhfoad/GCA_000195955.2_ASM19595v2_genomic.gff?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/GCA_000195955.2_ASM19595v2_genomic.gff
```

#### 1.6.2.2. Study Design File Download the study design file
The study design file is a tab separated file with 4 columns that encodes the phenotype and repeat information on the samples we will analyze in this practical. You will import this file into R to setup the differential expression analysis in the R Programming language. It is critical to keep track of the details of which samples are which etc. for downstream analysis purposes.  
```bash
wget https://www.dropbox.com/s/6y3z9btz3bg20pn/practical_study_design.txt?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/practical_study_design.txt
```

>**Note:** In this practical, we will use the existing study design file provided. However, in your own work, you will have to create this file on your own. You can create this file in several ways some of which could include exporting an Excel spreadsheet as a "Tab delimited" text file. You can also use your favorite text editor.such as using the **Text Editor** app under the Applications menu on the virtual machine.

<p align="center">
		<img src="https://github.com/WCSCourses/NGS_Bio_Africa/blob/main/images/H3ABioNet_Logo%20(1).png" style="width:100%">
</p>

#### 1.6.2.3. RData
Download `DE_data.RData` for later use. These files will be used for differential expression analysis in a later section of this tutorial.  
```bash
wget https://www.dropbox.com/s/6onagsnpp9wrkrv/DE_data.RData.zip?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/DE_data.RData.zip
```

Uncompress the `DE_data.RData.zip` file:
```bash
unzip DE_data.RData.zip
```

### 1.6.3 Setup R Analysis Environment
**Internet Access Required**   

In this section you will download additional `R` Libraries, packages of analysis tools to support differential expression analysis in `R`.

Start R by typing the following command:  
```bash
R
```

The Bioconductor Project provides tools written in R for the ’analysis and comprehension of high- throughput genomic data’. Here we will install two additional software packages (GenomicFeatures and tximport) for our practical today (See additional guidelines for Bioconductor package installation).

**Bioconductor Packages:**  
```r
install.packages("BiocManager")
BiocManager::install("GenomicFeatures", force = TRUE)
BiocManager::install("tximport")
```

**CRAN Packages:**  
```r
install.packages("pheatmap")
```

**When asked whether to:**  
`Update all/some/none? [a/s/n]` - select `no (n)`

**Exit `R`:**  
```r
q()
```

# 2. Introducing The Tutorial Dataset
Working through this tutorial, you will investigate the differences in expression between two Mycobacterium tuberculosis isolates with different phenotypes (**Sensitive** vs. **Resistant**). You will analyze 6 samples across the two conditions. The experiment is setup in the following way:

<div align="center">
	<table>
		<tr>
			<th> run </th>
			<th> unique_id </th>
			<th> phenotype </th>
			<th> repeat </th>
		</tr>
		<tr>
			<td> N2 </th>
			<td> sample1 </th>
			<td> Resistant </th>
			<td> 1 </th>
		</tr>
		<tr>
			<td> N6 </th>
			<td> sample2 </th>
			<td> Resistant </th>
			<td> 2 </th>
		</tr>
		<tr>
			<td> N10 </th>
			<td> sample3 </th>
			<td> Resistant </th>
			<td> 3 </th>
		</tr>
		<tr>
			<td> N14 </th>
			<td> sample4 </th>
			<td> Sensitive </th>
			<td> 1 </th>
		</tr>
		<tr>
			<td> N18 </th>
			<td> sample5 </th>
			<td> Sensitive </th>
			<td> 2 </th>
		</tr>
		<tr>
			<td> N22 </th>
			<td> sample6 </th>
			<td> Sensitive </th>
			<td> 3 </th>
		</tr>
    </table>
</div>

**Research Question:** what genes are differentially expressed between these two isolates that can explain the differences in their phenotypes?

**Check that you can see the FASTQ files in the practical directory:**  
```bash
ls N*.fq.gz
```

The FASTQ files contain the raw sequence reads for each sample. There are four lines per read: 
<ol>
	<li> Header
	<li> Sequence
	<li> Separator (usually a ‘+’)
	<li> Encoded quality value
</ol>

**Take a look at one of the FASTQ files:**  
```bash
zless N2_sub_R1.fq.gz | head
```

Find out more about FASTQ formats at <a href="https://en.wikipedia.org/wiki/FASTQ_format/">https://en.wikipedia.org/wiki/FASTQ_format/</a>.

## 2.1. Questions
**Q1: Why is there more than one FASTQ file per sample?** _Hint: think about why there is a `N2_sub_R1.fq.gz` and a `N2_sub_2.fq.gz`_

**Q2: How many reads were generated for the N2 sample?** _Hint: we want the total number of reads from both files (`N2_sub_R1.fq.gz` and `N2_sub_2.fq.gz`) so perhaps think about the FASTQ format and the number of lines for each read or whether there’s anything you can use in the FASTQ header to search and count._

**Q3: The three Resistant samples N2, N6 and N10 represent technical replicates. True or False? Comment on your answer.**

# 3. Estimate Transcript Abundance With `Salmon`
In this section, you will use <a href="https://salmon.readthedocs.io/en/latest/salmon.html">Salmon</a>, a transcript quantification method to estimate the number of reads in each sample that map to reference transcripts. This count is an estimate of the abundance or expression level of these transcripts in our experiment. As `Salmon` is an alignment free method, read quantification does not require an input BAM file. `Salmon` using a selective mapping algorithm to align the reads directly to a set of target transcripts such as those from a reference database for your organism.

Inputs include: 
<ul>
    <li> A set of target transcripts - FASTA format Reference Transcriptome GCA_000195955.2_ASM19595v2_genomic.transcripts.fa
    <li> Sample reads - FASTA/FASTQ files for your sample N2_sub_R1.fq.gz.....
</ul>

See more details in:
>Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. **Salmon provides fast and bias-aware quantification of transcript expression.** Nat Methods. 2017 Apr;14(4):417- 419. <a href="doi:10.1038/nmeth.4197">doi:10.1038/nmeth.4197</a>. PMID: 28263959; PMCID: PMC5600148.

## 3.1. Create Transcriptome Index
First, we create the necessary index files for any alignment tools downstream. The salmon index is created from the transcripts file, and while it is named `transcripts_index` here, you can name it anything.  
```bash
salmon index -t GCA_000195955.2_ASM19595v2_genomic.transcripts.fa -i transcripts_index -k 31
```

Take a look at the various files created in the `transcripts_index` folder created by `salmon`:

## 3.2. Transcript Quantification
Next we perform quantification with salmon using the transcript index folder we just created. However, it is important to perform key quality control analysis of your sample reads (`*.fq.gz`) files before proceeding with quantification.

### 3.2.1 Read Quality Control Analysis
Here we demonstrate how to trim and filter the reads using `Trimmomatic`. Here we clean the reads by looking at various quality metrics including:
<ul>
    <li> low quality reads - filter based on a min leading base quality of 3, and trailing base quality of 3. The scan with a sliding window that cuts when the average quality is lower than 15.
    <li> Illumina adapter sequences - clip out these technical artifacts using TruSeq3-PE.fa
    <li> read length - discard any ready with a length lower than 36
</ul>

```bash
trimmomatic PE N2_sub_R1.fq.gz N2_sub_R2.fq.gz \
    N2_trimmed_forward_paired.fq.gz N2_trimmed_forward_unpaired.fq.gz \
    N2_trimmed_reverse_paired.fq.gz N2_trimmed_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/manager/miniconda/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36
```

In a typically workflow, you would do `FastQC` analysis before and after to check the effects of the trimming, but as that is not the purpose of this section, we will skip that for now.

## 3.2.2 Transcript Quantificaiton using `salmon`
Now we used the trimmed and paired reads to estimate transcript abundance:  
```bash
salmon quant --geneMap GCA_000195955.2_ASM19595v2_genomic.gtf \
    --threads 2 -l A \
    -i transcripts_index GCA_000195955.2_ASM19595v2_genomic.fa \
    -1 N2_trimmed_forward_paired.fq.gz \
    -2 N2_trimmed_reverse_paired.fq.gz -o N2
```

This creates a new folder named “N2” in the directory, which contains a number of files:

<p align="center">
		<img src="https://github.com/WCSCourses/NGS_Bio_Africa/blob/main/images/H3ABioNet_Logo%20(1).png" style="width:100%">
</p>

The `quant` files are the files containing the read counts for the genes. These are what downstream tools like `DESeq2` or `edgeR` would use for the differential expression analysis.

It is not always practical to run commands for each sample one at a time. So we can write small scripts to process multiple samples using the same set of commands as is shown here:  
```bash
for r1 in *_R1.fq.gz
do
    echo $r1
    sample=$(basename $r1)
    sample=${sample%_sub_R1.fq.gz}
    echo "Processing sample: "$sample
    trimmomatic PE ${sample}_sub_R1.fq.gz ${sample}_sub_R2.fq.gz ${sample}_trimmed_forward_paired.fq.gz ${sample}_trimmed_forward_unpaired.fq.gz \
    ${sample}_trimmed_reverse_paired.fq.gz ${sample}_trimmed_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/home/manager/miniconda/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36

salmon quant --geneMap GCA_000195955.2_ASM19595v2_genomic.gtf \
    --threads 2 -l A -i transcripts_index GCA_000195955.2_ASM19595v2_genomic.fa \
    -1 ${sample}_trimmed_forward_paired.fq.gz \
    -2 ${sample}_trimmed_reverse_paired.fq.gz \
    -o ${sample}
done
```

## 3.3. Questions
In the `N2` folder, take a look at the `quant.sf` file and answer the following questions:

**Q1: What does TPM stand for?**

**Q2: How many reads in total are mapped to the tRNA genes?** _Hint: They start with “rna-”_

**Q3: Do you think this level of tRNA is acceptable?**

Take a moment to discuss what the effect of large amounts of reads from tRNA and rRNA could have on normalisation. _Hint: For further help understanding the quant.sf* output files look at this salmon documentation page._

# 4 Differential Expression Analysis with `DESeq2`
In this section, you will use the `R` Package `DESeq2` to perform differential expression analysis using the salmon quantification files you generated in the previous section. For more information about `DESeq` can look at the manuscript describing the method:
>Love MI, Huber W, Anders S (2014). **Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.** Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8.

We will now move to working in the `R` Programming environment. Start `R`:  
```bash
R
```

**Load the required libraries:**
>`R` has a base set of classes and methods and tools. The additional packages we installed are a set of functions and tools designed to handle specific biological data types and support analyses and visualization such as of RNA-seq data. Here, we load those software packages that are relevant to our analysis.

```r
library("tximportData")
library("tximport")
library("GenomicFeatures")
library("DESeq2")
library("pheatmap")
```

**Load the study design table:**   
```r
design_file <-"practical_study_design.txt"
samples <- read.table(design_file, header=TRUE)
```

**Set the row names as the samples names:**   
```r
rownames(samples) <- samples$run
samples
```
If we look at the samples object, we can see the data from our study design file.

**Link Salmon Output File Paths to Samples** - Create a files object pointing to the related `quant.sf` file for each sample:   
```r
root_dir <- "/home/manager/course_data/rna_seq_pathogen/practical"
files <- file.path(root_dir, samples$run, "quant.sf")
files
```

In the `files` object, we should now see the path to the `quant.sf` file in each sample’s folder. This is much simpler and less error-prone than typing out each file path manually.

**Use the run column to map samples to their file paths:**   
```r
names(files) <- samples$run
```

The files object now has the sample name linked to the file path and can be used as a way to map the information.

### 4.0.1 Make Annotation Database Object

**Load the annotation information:**   
```r
path_to_gff <- "GCA_000195955.2_ASM19595v2_genomic.gff"
txdb <-makeTxDbFromGFF(path_to_gff, organism="Mycobacterium tuberculosis")
```

The `makeTxDbFromGFF` function is part of the `GenomicFeatures` library and is used to create a Transcript Database (`txdb`) object from a GFF annotation file. The `TxDb` class is a container for storing transcript annotations.

**Extract the GENEID key values from the tbdx object:**   
```r
k <- keys(txdb, keytype = "GENEID")
head(k)
```

The `k` object is a list of all the gene names based on the GENEID as extracted from the annotation file.

**Create a mapping table based on the GENEID and TXNAME info:**   
```r
tx2gene <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
head(tx2gene)
```

Take a look at the output, you will see 2 columns. This will be used to map the gene names from the salmon files to the annotation files.

**Rename the genes in the GENEID column** - We do this to match the gene names found in the salmon `quant.sf` files with those found in the transcriptome and annotation files:  
```r
tx2gene[["GENEID"]] <- with(tx2gene, ifelse(!grepl("rna-", TXNAME),paste0("gene-", TXNAME), TXNAME))
head(tx2gene[["GENEID"]])
```

The above command goes through the `tx2gene` object, and looks for where the gene name in the **TXNAME** column DOES NOT have the prefix "rna-". In those rows, it adds the "gene-" prefix to the gene name in the **GENEID** column.

When the annotation file is imported by `makeTxDbFromGFF`, it removes the "gene-" prefix that is present in the `quant` files. This is because of the way that the annotation and transcript files are formatted.

You will see in the GFF and transcript files, that the gene ID is formatted like this: `ID=gene-Rv0005`  
Whereas in the GTF file it is formatted like this: `gene_id ”Rv0005”`

Though salmon uses the GTF file, it adds on the gene- prefix where it is missing while `makeTxDbFromGFF` removes it where it is present.  
```r
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
```

### 4.0.2 Create the DESeq data set object
```r
ddsTxi <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ phenotype)
```

The inputs to create a DESeq object are:
<ul>
    <li> txi - the summary of transcript level abundance estimates
    <li> colData - information from the study design file loaded earlier
    <li> design - formula that expresses the relationship between the gene counts and the variables in the study design. Here, we are using phenotype.
</ul>

The design formula is a statement that specifies how we want to model the variation in gene expression from the abundance estimates for each sample (counts). This statement is loosely saying that we expect that the level of expression of a gene is dependent on the phenotype of the bacterial isolate. If we have more than a single experimental factor to consider, we would change how we specify the design formula.

>**Note:** For a more extensive treatment of how to setup the design formula for more complex experi- mental designs, read through A guide to creating design matrices for gene expression experiments.

### 4.0.3 Perform DESeq Analysis
**Normalisation and model fitting with DESeq2:**  
Next we used `DESeq2` to normalize the data and fit a model that relates the gene count information to the phenotype.  
```r
ddsTxi <- DESeq(ddsTxi)
```

This step does the actual normalisation and model fitting for the data, and results in a `deseq` data set.

**Load provided `.RData` for differential expression analysis:**   
```r
load("DE_data.RData")
```

## 4.1. Exploratory Visualization
**Transform data for visualization:**   
The `rlog` transformation provides functionality for visualization and clustering.  
```r
rldTxi <- rlog(ddsTxi, blind = FALSE)
```

**Plot principal component analysis (PCA):**  
The PCA plot allow us to assess the similarities and differences between the study samples in order to determine if the data fit our expectation compared to the experimental design. This is a good quality control check to use to ensure that we did not introduce any errors during sample processing, sequencing and primary data analysis steps.  
```r
plotPCA(rldTxi, intgroup = c("phenotype"))
```

## 4.2. Result Generation and Exploration
**Get summary of differentially expression results:**  
The results function of `DESeq2` is used to get the the results of the `DESeq` function ranked by metrics to use to determine which genes are significantly differentially expressed.  
```r
summary(results(ddsTxi))
```

**Apply LFC threshold and adjusted p-value cutoffs to get significantly differentially expressed genes:**   
A significance threshold (`FDR`, or `alpha`) and `log2FC` threshold can be passed to the results function to filter non-signficant differentially expressed genes.  
```r
resTxi <- results(ddsTxi, alpha=0.05, lfcThreshold=0.5,altHypothesis="greaterAbs")
summary(resTxi)
```

**Get significantly differentially expressed genes:**   
```r
sigTxi <- resTxi[!is.na(resTxi$padj) & resTxi$padj < 0.05 & abs(resTxi$log2FoldChange) >= 0.5,]
sigTxi
```

**Visualize your differentially expressed genes:**  
Extract the `rlog` transformed expression values for the significant gene to use for the heatmap.  
```r
matrixTxi <- assay(rldTxi)[rownames(sigTxi), ]
matrixTxi
```

**Plot heatmap:**  
```r
pheatmap(matrixTxi, cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=TRUE)
```

## 4.3. Questions
**Q1. How many genes are up- and down-regulated between the resistant and sensitive isolates of the Mycobacterium tuberculosis?** _Hint:use the summary function on the DESeq results object._

**Q2: How many genes are significantly differentially expressed (i.e., meet the LFC threshold and adjusted p-value cutoffs) between the resistant and sensitive isolates of the Mycobacterium tuberculosis? Name these genes.**

**Q3. What are the _p_-values for the significantly differentially expressed genes?**

# 5. Key Aspects of Differential Expression Analysis
## 5.1. Replicates and power
**Biological replicates** are parallel measurements of biologically distinct samples that capture random biological variation, which may itself be a subject of study or a noise source.

**Technical replicates** are repeated measurements of the same sample that represent independent measures of the random noise associated with protocols or equipment

>Blainey, Paul et al. **“Points of significance: replication.”** Nature methods vol. 11,9 (2014): 879-80. <a href="doi:10.1038/nmeth.3091">doi:10.1038/nmeth.3091</a>

In order to accurately ascertain which genes are differentially expressed and by how much it is necessary to use replicated data. As with all biological experiments doing it once is simply not enough. There is no simple way to decide how many replicates to do, it is usually a compromise between statistical power and cost. By determining how much variability there is in the sample preparation and sequencing reactions, we can better assess how highly genes are really expressed and more accurately determine any differences. The key to this is performing biological rather than technical replicates. This means, for instance, growing up three batches of parasites, treating them all identically, extracting RNA from each and sequencing the three samples separately. Technical replicates, whereby the same sample is sequenced three times do not account for the variability that really exists in biological systems or the experimental error between batches of parasites and RNA extractions.

>**Note:** More replicates will help improve power for genes that are already detected at high levels, while deeper sequencing will improve power to detect differential expression for genes which are expressed at low levels.

## 5.2. _p_-values vs. _q_-values
When asking whether a gene is differentially expressed we use statistical tests to assign a _p_-value. If a gene has a _p_-value of 0.05, we say that there is only a 5% chance that it is not really differentially expressed. However, if we are asking this question for every gene in the genome, then we would expect to see _p_-values less than 0.05 for many genes even though they are not really differentially expressed. Due to this statistical problem, we must correct the _p_-values so that we are not tricked into accepting a large number of erroneous results. _q_-values are _p_-values which have been corrected for what is known as multiple hypothesis testing. Therefore, it is a qvalue of less than 0.05 that we should be looking for when asking whether a gene is signficantly differentially expressed.

## 5.3. Altermate software
If you have a good quality genome and genome annotation such as for model organisms e.g. human, mouse, _Plasmodium_ etc. map to the transcriptome to determine transcript abundance. This is even more relevant if you have variant transcripts per gene as you need a tool which will do its best to determine which transcript is really expressed. Kallisto (<a href="https://pachterlab.github.io/kallisto/">Bray et al. 2016; PMID: 27043002</a>) and eXpress (<a href="https://pubmed.ncbi.nlm.nih.gov/23160280/">Roberts & Pachter, 2012; PMID: 23160280</a>) are examples of these tools.

## 5.4. What do I do with a list of differentially expressed genes?
Differential expression analysis results are a list of genes which show differences between conditions. It can be daunting trying to determine what the results mean. On one hand, you may find that there are no real differences in your experiment. Is this due to biological reality or noisy data? On the other hand, you may find several thousands of genes are differentially expressed.

> What can you say about that?

Other than looking for genes you expect to be different or unchanged, one of the first things to do is look at common themes in gene functionality across your list. For example, you can carry out Gene Ontology (GO) term enrichment analysis. There are many different algorithms for this, but you could annotate your genes with functional terms from GO using for instance Blast2GO (<a href="https://pubmed.ncbi.nlm.nih.gov/16081474/">Conesa et al., 2005; PMID: 16081474</a>) and then use TopGO (<a href="https://pubmed.ncbi.nlm.nih.gov/16606683/">Alexa et al., 2005; PMID: 16606683</a>) to determine whether any particular sorts of genes occur more than expected in your differentially expressed genes.

<h1 align="center">Congratulations, you have reached the end of this tutorial!</h1>

# 6. Normalisation
## 6.1 Introduction
In the previous section, we looked at estimating transcript abundance with `Kallisto`. The abundances are reported as transcripts per million (TPM), but what does TPM mean and how is it calculated?

The objectives of this part of the tutorial are:
<ul>
	<li> understand why RNA-Seq normalisation metrics are used
	<li> understand the difference between RPKM, FPKM and TPM
	<li> calculate RPKM and TPM for a gene of interest
</ul>

There are many useful websites, publications and blog posts which go into much more detail about RNA-Seq normalisation methods. Here are just a couple (in no particular order):
<ul>
	<li> <a href="https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/">What the FPKM? A review of RNA-Seq expression units</a>
	<li> <a href="https://statquest.org/2015/07/09/rpkm-fpkm-and-tpm-clearly-explained/">RPKM, FPKM and TPM, clearly explained</a>
	<li> <a href="https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4728800/">A survey of best practices for RNA-seq data analysis</a>
	<li> <a href="http://robpatro.com/blog/?p=235">The RNA-seq abundance zoo</a>
</ul>

## 6.2. Why do we use normalisation units instead of raw counts?
Raw reads counts are the number of reads originating from each transcript which can be affected by several factors:
<ul>
	<li> <b>Sequencing depth (total number of reads)</b><br>The more we sequence a sample, the more reads we expect to be assigned.
	<li> <b>Gene/transcript length</b><br>The longer the gene or transcript, the more reads we expect to be assigned to it.
</ul>

<p align="center">
    <img src="https://github.com/WCSCourses/NGS_Bio_Africa/blob/main/images/H3ABioNet_Logo%20(1).png" style="width:100%">
    <b>Figure 4.</b> Effect of sequencing depth and gene length on raw read counts.
</p>

>Look at the top part of Figure 4. In which sample, X or Y, is the gene more highly expressed?

Neither, it’s the same in both. What we didn’t tell you was that the total number of reads generated for sample A was twice the number than for sample B. That meant almost twice the number of reads are assigned to the same gene in sample A than in sample B.

>Look at the bottom part of Figure 4. Which gene, X or Y, has the greatest gene level expression?

Neither, they are both expressed at the same level. This time we didn’t tell you that gene X is twice the length of gene Y. This meant that almost twice the number reads were assigned to gene X than gene Y.

In the top part of **Figure 4**, the gene in sample X has twice the number of reads assigned to it than the same gene in sample Y. What isn’t shown is that sample X had twice the number or total reads than sample Y so we would expect more reads to be assigned in sample X. Thus, the gene is expressed at roughly the same level in both samples. In the bottom part of **Figure 4**, gene X has twice the number of reads assigned to it than gene Y. However, gene X is twice the length of gene Y and so we expect more reads to be assigned to gene X. Again, the expression level is roughly the same.

### 6.2.1. Reads per kilobase per million (RPKM)
Reads per kilobase (of exon) per million (reads mapped) or RPKM is a within sample normalisation method which takes into account sequencing depth and length biases.

To calculate RPKM, you first normalise by sequencing depth and then by gene/transcript length.
<ol>
    <li> <b>Get your per million scaling factor</b> <br> Count up the total number of reads which have been assigned (mapped) in the sample. Divide this number by 1,000,000 (1 million) to get your per million scaling factor (N).
	<li> <b>Normalise for sequencing depth</b> <br> Divide the number of reads which have been assigned to the gene or transcript (C) by the per million scaling factor you calculated in step 1. This will give you your reads per million (RPM).
	<li> <b>Get your per kilobase scaling factor</b> <br> Divide the total length of the exons in your transcript or gene in base pairs by 1,000 (1 thousand) to get your per kilobase scaling factor (L).
	<li> <b>Normalise for length</b> <br> Divide your RPM value from step 2 by your per kilobase scaling factor (length of the gene/transcript in kilobases) from step 3. This will give you your reads per kilobase per million or RPKM.
</ol>

This can be simplified into the following equation:

<p align="center">
	<img src="https://render.githubusercontent.com/render/math?math=RPKM = \frac{C}{LN}" style="width:25%">
</p>

Where:
<ul>
	<li> C is number of reads mapped to the transcript or gene
	<li> L is the total exon length of the transcript or gene in kilobases
	<li> N is the total number of reads mapped in millions
</ul>

### 6.2.2 Fragments per kilobase per million (FPKM)
Fragments per kilobase per million or FPKM is essentially the same as RPKM except that: 
<ul>
    <li> RPKM is designed for single-end RNA-Seq experiments
    <li> FPKM is designed for paired-end RNA-Seq experiments
</ul>

In a paired-end RNA-Seq experiment, two reads may be assigned to a single fragment (in any orientation). Also, in some cases, only one of those reads will be assigned to a fragment (singleton). The only difference between RPKM and FPKM is that FPKM takes into consideration that two reads may be assigned to the same fragment.

### 6.2.3 Transcripts per million (TPM)
Calculating the transcripts per million or TPM is a similar process to RPKM and FPKM. The main difference is that you will first normalise for length bias and then for sequencing depth bias. In a nut shell, we are swapping the order of normalisations.
<ol>
	<li> <b>Get your per kilobase scaling factor</b> <br> Divide the total length of the exons in your transcript in base pairs by 1,000 (1 thousand) to get your per kilobase scaling factor.
	<li> <b>Normalise for length</b> <br> Divide the number of reads which have been assigned to the transcript by the per kilobase scaling factor you calculated in step 1. This will give you your reads per kilobase (RPK).
	<li> <b>Get the sum of all RPK values in your sample</b> <br> Calculate the RPK value for all of the transcripts in your sample. Add all of these together to get your total RPK value.
	<li> <b>Get your per million scaling factor</b> <br> Divide your total RPK value from step 3 by 1,000,000 (1 million) to get your per million scaling factor.
	<li> <b>Normalise for sequencing depth</b> <br> Divide your RPK value calculated in step 2 by the per million scaling factor from step 4. You now have your transcripts per millions value or TPM.
</ol>

## 6.3. Calculating RPKM and TPM values
To try and answer this, let’s look at a worked example. Here, we have three genes (A-C) and three biological replicates (1-3).

<div align="center">
    <table>
        <tr>
            <th> Gene </th>
			<th> Length </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 2000 bases </td>
			<td> 10 </td>
            <td> 12</td>
            <td> 30 </td>
		<tr>
			<td> B </td>
			<td> 4000 bases </td>
			<td> 20 </td>
            <td> 25 </td>
            <td> 60 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1000 bases </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

There are two things to notice in our dataset:
<ul>
	<li> Gene B has twice number reads mapped than gene A, possibly as it’s twice the length
	<li> Replicate 3 has more reads mapped than any of the other replicates, regardless of which gene we look at.
</ul>

### 6.3.1. Calculating RPKM
#### Step 1: get your per million scaling factor

In the table below is the total number of reads which mapped for each of the replicates. To get our per million scaling factor, we divide each of these values by 1,000,000 (1 million).

<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
			<th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> Total reads mapped  </td>
			<td> 3 500 000 </td>
			<td> 4 500 000 </td>
            <td> 10 600 000</td>
		<tr>
			<td> Per million reads </td>
			<td> 3.5 </td>
			<td> 4.5 </td>
            <td> 10.6 </td>
        </tr>
    </table>
</div>

#### Step 2: normalise for sequencing depth
We now divide our read counts by the per million scaling factor to get our reads per million (RPM).

**Before:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 10 </td>
            <td> 12</td>
            <td> 30 </td>
		<tr>
			<td> B </td>
			<td> 20 </td>
            <td> 25 </td>
            <td> 60 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

**After:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 (RPM) </th>
            <th> Replicate 2 (RPM) </th>
            <th> Replicate 3 (RPM) </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 2.857 </td>
            <td> 2.667 </td>
            <td> 2.830 </td>
		<tr>
			<td> B </td>
			<td> 5.714 </td>
            <td> 5.556 </td>
            <td> 5.660 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1.429 </td>
            <td> 1.778 </td>
            <td> 1.415 </td>
        </tr>
    </table>
</div>

#### Step 3: get your per kilobase scaling factor
Here we have our gene length in base pairs. For our per kilobase scaling factor we need to get our gene length in kilobases by dividing it by 1,000.

<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Length (base pairs) </th>
            <th> Length (kilobases) </th>
		</tr>
		<tr>
			<td> A </td>
            <td> 2,000 </td>
            <td> 2 </td>
		<tr>
			<td> B </td>
			<td> 4,000 </td>
            <td> 4 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1,000 </td>
            <td> 1 </td>
        </tr>
    </table>
</div>

#### Step 4: normalise for length
Finally, we divide our RPM values from step 2 by our per kilobase scaling factor from step 3 to get our reads per kilobase per million (RPKM).

**Before:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 (RPM) </th>
            <th> Replicate 2 (RPM) </th>
            <th> Replicate 3 (RPM) </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 2.857 </td>
            <td> 2.667 </td>
            <td> 2.830 </td>
		<tr>
			<td> B </td>
			<td> 5.714 </td>
            <td> 5.556 </td>
            <td> 5.660 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1.429 </td>
            <td> 1.778 </td>
            <td> 1.415 </td>
        </tr>
    </table>
</div>

**After:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 (RKPM) </th>
            <th> Replicate 2 (RKPM) </th>
            <th> Replicate 3 (RKPM) </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 1.43 </td>
            <td> 1.33 </td>
            <td> 1.42 </td>
		<tr>
			<td> B </td>
			<td> 1.43 </td>
            <td> 1.39 </td>
            <td> 1.42 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1.43 </td>
            <td> 1.78 </td>
            <td> 1.42 </td>
        </tr>
    </table>
</div>

Notice that even though replicate 3 had more reads assigned than the other samples and a greater sequencing depth, its RPKM is quite similar. And, that although gene B had twice the number of reads assigned than gene A, its RPKM is the same. This is because we have normalised by both length and sequencing depth.

### 6.3.2. Calculating TPM
Now we’re going to calculate the TPM values for the same example data. As a reminder, here are our three genes (A-C) and three biological replicates (1-3).

<div align="center">
    <table>
        <tr>
            <th> Gene </th>
			<th> Length </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 2000 bases </td>
			<td> 10 </td>
            <td> 12</td>
            <td> 30 </td>
		<tr>
			<td> B </td>
			<td> 4000 bases </td>
			<td> 20 </td>
            <td> 25 </td>
            <td> 60 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1000 bases </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

#### Step 1: get your per kilobase scaling factor
Again, our gene lengths are in base pairs. For our per kilobase scaling factor we need to get our gene length in kilobases by dividing it by 1,000.

<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Length (base pairs) </th>
            <th> Length (kilobases) </th>
		</tr>
		<tr>
			<td> A </td>
            <td> 2,000 </td>
            <td> 2 </td>
		<tr>
			<td> B </td>
			<td> 4,000 </td>
            <td> 4 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1,000 </td>
            <td> 1 </td>
        </tr>
    </table>
</div>


#### Step 2: normalise for length
Now we divide the number of reads which have been assigned to each gene by the per kilobase scaling factor we just calculated. This will give us our reads per kilobase (RPK).

**Before:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 10 </td>
            <td> 12</td>
            <td> 30 </td>
		<tr>
			<td> B </td>
			<td> 20 </td>
            <td> 25 </td>
            <td> 60 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

**After:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 5 </td>
            <td> 6 </td>
            <td> 15 </td>
		<tr>
			<td> B </td>
			<td> 5 </td>
            <td> 6.25 </td>
            <td> 15 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

#### Step 3: get the sum of all RPK values in your sample
Next, we sum the RPK values for each of our replices. This will give use our total RPK value for each replicate. To make this example scalable, we assume there are other genes so the total RPK is made up.

<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 5 </td>
            <td> 6 </td>
            <td> 15 </td>
		<tr>
			<td> B </td>
			<td> 5 </td>
            <td> 6.25 </td>
            <td> 15 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
        <tr>
	        <td> ... </td>
			<td> ... </td>
			<td> ... </td>
            <td> ... </td>
        </tr>
		<tr>
			<th> Total RPK </th>
			<th> 150 000 </th>
			<th> 202 500 </th>
            <th> 450 000 </th>
        </tr>
    </table>
</div>

#### Step 4: get your per million scaling factor
Here, instead of dividing our total mapped reads by 1,000,000 (1 million) to get our per million scaling factor, we divide our total RPK values by 1,000,000 (1 million).

<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> Total RPK </td>
			<td> 150 000 </td>
            <td> 202 500 </td>
            <td> 450 000 </td>
		<tr>
			<td> Per million RPK </td>
			<td> 0.1500 </td>
            <td> 0.2025 </td>
            <td> 0.4500 </td>
        </tr>
    </table>
</div>

#### Step 5: normalise for sequencing depth
Finally, we divide our individual RPK values from step 2 by the per million scaling factor in step 4 to give us our TPM values. 

**Before:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 5 </td>
            <td> 6 </td>
            <td> 15 </td>
		<tr>
			<td> B </td>
			<td> 5 </td>
            <td> 6.25 </td>
            <td> 15 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 5 </td>
            <td> 8 </td>
            <td> 15 </td>
        </tr>
    </table>
</div>

**After:**
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 33.33 </td>
            <td> 29.63 </td>
            <td> 33.33 </td>
		<tr>
			<td> B </td>
			<td> 33.33 </td>
            <td> 30.86 </td>
            <td> 33.33 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 33.33 </td>
            <td> 39.51 </td>
            <td> 33.33 </td>
        </tr>
    </table>
</div>

## 6.4. Which normalisation unit should I use?
Well, there’s a lot of debate around this, so let’s look at our total normalised values for each replicate.

### 6.4.1. RPKM
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 (RKPM) </th>
            <th> Replicate 2 (RKPM) </th>
            <th> Replicate 3 (RKPM) </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 1.43 </td>
            <td> 1.33 </td>
            <td> 1.42 </td>
		<tr>
			<td> B </td>
			<td> 1.43 </td>
            <td> 1.39 </td>
            <td> 1.42 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 1.43 </td>
            <td> 1.78 </td>
            <td> 1.42 </td>
        </tr>
        <tr>
            <th> Total RPKM </th>
            <th> 4.29 </th>
            <th> 4.50 </th>
            <th> 4.25 </th>
        </tr>
    </table>
</div>

### 6.4.2. TPM
<div align="left">
    <table>
        <tr>
            <th> Gene </th>
			<th> Replicate 1 </th>
            <th> Replicate 2 </th>
            <th> Replicate 3 </th>
		</tr>
		<tr>
			<td> A </td>
			<td> 33.33 </td>
            <td> 29.63 </td>
            <td> 33.33 </td>
		<tr>
			<td> B </td>
			<td> 33.33 </td>
            <td> 30.86 </td>
            <td> 33.33 </td>
        </tr>
		<tr>
			<td> C </td>
			<td> 33.33 </td>
            <td> 39.51 </td>
            <td> 33.33 </td>
        </tr>
        <tr>
            <th> Total TPM </th>
            <th> 100 </th>
            <th> 100 </th>
            <th> 100 </th>
        </tr>
    </table>
</div>

Notice that that total TPM value for each of the replicates is the same. This is not true for RPKM and FPKM where the total values differ. With TPM, having the same total value for each replicate makes it easier to compare the proportion of reads mapping to each gene across replicates (although you shouldn’t really compare across experiments). With RPKM and FPKM, the differing total values make it much harder to compare replicates.
