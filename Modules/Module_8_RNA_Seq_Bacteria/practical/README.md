# RNA-Seq Expression Analysis

## 1.1 Introduction
RNA sequencing (RNA-Seq) is a high-throughput method used to profile the transcriptome, quantify gene expression and discover novel RNA molecules. This tutorial uses RNA sequencing of from two _**Mycobacterium tuberculosis**_ isolates (each with 3 replicates - 6 total samples) to walk you through transcriptome pseudo-alignment and quantification and how to profile transcriptomic differences between experimental conditions (**Sensitive** vs. **Resistant**) by identifying differentially expressed genes.

For an introduction to RNA-Seq principles and best practices see:

> **A survey of best practices for RNA-Seq data analysis.** Ana Conesa, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, Daniel J. Gaffney, Laura L. Elo, Xuegong Zhang and Ali Mortazavi Genome Biol. 2016 Jan 26;17:13 doi:10.1186/s13059-016-0881-8

## 1.2. Learning Outcomes

By the end of this practical, you should be familiar with:
<ul>
    <li> Using salmon to quantify transcript abundance 
    <li> Knowing how to create a study design file
    <li> Loading the required files into R
    <li> Creating a DESeq2 data set
    <li> Running DESeq
    <li> Visualising the relationship between samples 
    <li> Creating a results object
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
    <li> Quantification - How many reads come from each gene? (<b>Salmon</b>) 
    <li> Normalisation - Dealing with biases in the data (<b>DESeq2</b>)
    <li> Differential expression analysis (<b>DESeq2</b>)
    <li> Visualisation and reporting (<b>R libraries</b>)
</ul>

**Additional Resources**   
<a href="https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html">RNA-seq workflow: gene-level exploratory analysis and differential expression</a>

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

```
cd /home/manager/course_data/rna_seq_pathogen
```

Create a working directory for the practical:

```
mkdir practical
```

Create a copy of the practial dataset to maintain the data integrity of the original dataset:

```
 cp /home/manager/course_data/rna_seq_pathogen/data/bacterial/* /home/manager/course_data/rna_seq_pathogen/practical
```

Copy and paste the command above as a single line at your command line window.   
>Note: if you make a mistake at any point during this tutorial, you can reset by deleting the practical folder and restarting from this section.

Move to the practical directory:

```
cd practical
```

### 1.6.2. Download Supplemental Datasets

**Internet Access Required**

#### Reference Transcriptome File

Download the GFF formatted version of the annotation file. Unfortunately due to the lack of standardization in bioinformatics, we need both the <a href="https://mblab.wustl.edu/GTF22.html">>GTF</a> file and the <a href="https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md">GFF</a> file for analysis (See the differences and similarities between these two file format specifications <a href="https://www.ensembl.org/info/website/upload/gff.html?redirect=no">here</a>). Salmon uses the GTF file and the GenomicFeatures library in R needs a GFF file. Both files contain the same annotation information, they are just formatted differently. No tools are able to reliably convert one to the other because even GFF files are not formatted consistently.

```
wget https://www.dropbox.com/s/4yjgbmy3dyhfoad/GCA_000195955.2_ASM19595v2_genomic.gff?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/GCA_000195955.2_ASM19595v2_genomic.gff
```

#### Study Design File Download the study design file

The study design file is a tab separated file with 4 columns that encodes the phenotype and repeat information on the samples we will analyze in this practical. You will import this file into R to setup the differential expression analysis in the R Programming language. It is critical to keep track of the details of which samples are which etc. for downstream analysis purposes.

```
wget https://www.dropbox.com/s/6y3z9btz3bg20pn/practical_study_design.txt?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/practical_study_design.txt
```

>**Note:** In this practical, we will use the existing study design file provided. However, in your own work, you will have to create this file on your own. You can create this file in several ways some of which could include exporting an Excel spreadsheet as a "Tab delimited" text file. You can also use your favorite text editor.such as using the **Text Editor** app under the Applications menu on the virtual machine.

<p align="center">
		<img src="https://github.com/WCSCourses/NGS_Bio_Africa/blob/main/images/H3ABioNet_Logo%20(1).png" style="width:100%">
</p>

#### RData

Download `DE_data.RData` for later use. These files will be used for differential expression analysis in a later section of this tutorial.

```
wget https://www.dropbox.com/s/6onagsnpp9wrkrv/DE_data.RData.zip?dl=1 -O /home/manager/course_data/rna_seq_pathogen/practical/DE_data.RData.zip
```

Uncompress the `DE_data.RData.zip` file:

```
unzip DE_data.RData.zip
```

### 1.6.3 Setup R Analysis Environment

**Internet Access Required**   
In this section you will download additional `R` Libraries, packages of analysis tools to support differential expression analysis in `R`.   
Start R by typing the following command:

```
R
```

The Bioconductor Project provides tools written in R for the ’analysis and comprehension of high- throughput genomic data’. Here we will install two additional software packages (GenomicFeatures and tximport) for our practical today (See additional guidelines for Bioconductor package installation).

**Bioconductor Packages:**

```
install.packages("BiocManager")
BiocManager::install("GenomicFeatures", force = TRUE)
BiocManager::install("tximport")
```

**CRAN Packages:**

```
install.packages("pheatmap")
```

**When asked whether to:**   
Update all/some/none? [a/s/n] - select no (n)

**Exit `R`:**

```
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

Check that you can see the FASTQ files in the practical directory:   
```
ls N*.fq.gz
```



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
