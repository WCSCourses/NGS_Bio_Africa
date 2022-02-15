# 1. RNA-Seq Expression Analysis

## 1.1. Introduction

RNA sequencing (RNA-Seq) is a high-throughput method used to profile the transcriptome, quan- tify gene expression and discover novel RNA molecules or splice variants. This tutorial uses RNA sequencing of human cancer samples to walk you through transcriptome alignment, visualisation, simple quality control checks and shows you how to profile transcriptomic differences by identifying differentially expressed genes.

For an introduction to RNA-Seq principles and best practices see:

> **A survey of best practices for RNA-Seq data analysis.** Ana Conesa, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, Daniel J. Gaffney, Laura L. Elo, Xue- gong Zhang and Ali Mortazavi Genome Biol. 2016 Jan 26;17:13 doi:10.1186/s13059-016-0881-8

## 1.2. Learning Outcomes

By the end of this tutorial you can expect to be able to:
<ul>
	<li> Align RNA-Seq reads to a reference genome and a transcriptome
	<li> Visualise transcription data using standard tools
	<li> Perform QC of NGS transcriptomic data
	<li> Quantify the expression values of your transcripts using standard tools
</ul>

## 1.3. Practical Outline

This tutorial comprises the following sections:  
<ol>
	<li>Introducing the tutorial dataset  
	<li>Mapping RNA-Seq reads to the genome with HISAT2  
	<li>Visualising transcriptomes with IGV  
	<li>Transcript quantification with Kallisto  
	<li>Identifying differentially expressed genes with Sleuth  
	<li>Key aspects of differential expression analysis
</ol>

## 1.4. Authors

This tutorial was developed by Victoria Offord and Adam Reid and adapted for use here by Nyasha Chambwe.

## 1.5. Prerequisites

This tutorial assumes that you have the following software or packages and their dependencies installed on your computer. The software or packages used in this tutorial may be updated from time to time, so we have given you the version which was used when writing the tutorial.

<div align="center">
	<table>
		<tr>
			<th width="20%"> Package </th>
			<th width="60%"> Link for download/installation instructions </th>
			<th width="20%"> Version tested</th>
		</tr>
		<tr>
			<td> HISAT2 </td>
			<td> <a href="https://ccb.jhu.edu/software/hisat2/index.shtml">https://ccb.jhu.edu/software/hisat2/index.shtml</a> </td>
			<td> 2.1.0 </td>
		</tr>
		<tr>
			<td> samtools </td>
			<td> <a href="https://github.com/samtools/samtools">https://github.com/samtools/samtools</a> </td>
			<td> 1.10 </td>
		</tr>
		<tr>
			<td> IGV </td>
			<td> <a href="https://software.broadinstitute.org/software/igv/">https://software.broadinstitute.org/software/igv/</a> </td>
			<td> 2.7.2 </td>
		</tr>
		<tr>
			<td> kallisto </td>
			<td> <a href="https://pachterlab.github.io/kallisto/download">https://pachterlab.github.io/kallisto/download</a> </td>
			<td> 0.46.2 </td>
		</tr>
		<tr>
			<td> R </td>
			<td> <a href="https://www.r-project.org/">https://www.r-project.org/</a> </td>
			<td> 4.0.2 </td>
		</tr>
		<tr>
			<td> sleuth </td>
			<td> <a href="https://pachterlab.github.io/sleuth/download">https://pachterlab.github.io/sleuth/download</a> </td>
			<td> 0.30.0 </td>
		</tr>
		<tr>
			<td> bedtools </td>
			<td> <a href="http://bedtools.readthedocs.io/en/latest/content/installation.html">http://bedtools.readthedocs.io/en/latest/content/installation.html</a> </td>
			<td> 2.29.2</td>
		</tr>
	</table>
</div>

## 1.6. Where can I find the tutorial data?

You can find the data for this tutorial by typing the following command in a new terminal window:

```
cd /home/manager/course_data/rna_seq_human
```

Within the module directory, create a directory to write the outputs you will generate during this practical using the following command:

```
mkdir -p outputs
```

### 1.6.1. Download Suppplemental Datasets

**Internet Access Required**  
We will need to download two additional annotation files for this tutorial. Follow these instructions:

**Use wget to download supplemental data files:**  

```
wget https://www.dropbox.com/s/nfuea7ik6wlidum/hsapiens_chr21_transcript_to_gene.csv?dl=1 -O data/hsapiens_chr21_transcript_to_gene.csv
wget https://www.dropbox.com/s/8xt8q1o0aej1ry1/hsapiens_chr21_transcripts.fa?dl=1 -O data/hsapiens_chr21_transcripts.fa
```

Confirm that you have the following files in the data folder:
<ul>
	<li> `hsapiens_chr21_transcript_to_gene.csv`
	<li> `hsapiens_chr21_transcripts.fa`
</ul>

```
ls -lhr data/hsapiens_chr21_transcript*
```

You are ready to go. Now, let’s head to the first section of this tutorial which will be introducing the tutorial dataset.

# 2. Introducing The Tutorial Dataset

Prostate cancer incidence and the risk of dying from the disease is higher in **African American (AA)** men compared to **European American (EA)** men (<a href="https://cancerprogressreport.aacr.org/disparities/">American Cancer Association for Cancer Research AACR Cancer Disparities Report 2020</a>). These disparities are thought to be caused by a complex interplay between environmental and lifestyle factors as well as biological factors related to ancestry genetics. The molecular basis for these disparities is not well understood.  

Working through this tutorial, you will analyse RNA-sequencing data from high grade prostate can- cer and matched noncancer adjacent tissue samples from a cohort of ethnically diverse men (AAs and EAs) to investigate the molecular differences between these groups.
The dataset you will be using for this tutorial have been taken from the following publication:
> **Exogenous IL-6 induces mRNA splice variant MBD2_v2 to promote stemness in TP53 wild-type, African American PCa cells.** Teslow, E. A., Bao, B., Dyson, G., Leg- endre, C., Mitrea, C., Sakr, W., Carpten, J. D., Powell, I., & Bollig-Fischer, A. (2018). Molecular oncology, 12(7), 1138–1152. <a href="https://doi.org/10.1002/1878-0261.12316">https://doi.org/10.1002/1878-0261.12316</a>

This published dataset is publicly available from the Gene Expression Omnibus Database (Accession Number: <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104131">GSE104131</a>) and has been modified and adapted for use in these education materials.

## 2.1. Research Question

This dataset can be used to address two research questions namely:
<ul>
	<li> To what extent can we determine race-specific differential expression in prostate cancer?
	<li>What genes are differentially expressed between tumour and normal prostate cancer samples, taking into account any race-specific effects?
</ul>

The hypothesis underlying this analysis is that ethnicity-based differences in disease outcomes are driven by genetics define ancestry groups and can be observed at the molecular level.

## 2.2. Excercise 1

To illustrate the principles of RNA-seq analysis, you will analyse 12 RNA-seq samples that represent tumour-normal pairs (**prostate tumor PT; normal normal prostate NP**) from six individuals. Three individuals are African American (AA) and three are European American (EA).

<div align="center">
	<table>
		<tr>
			<th> patientID </th>
			<th> normalID </th>
			<th> tumorID </th>
			<th> ethnicity </th>
		</tr>
		<tr>
			<td> 487 </td>
			<td> NP4 </td>
			<td> PT4 </td>
			<td> AA </td>
		</tr>
		<tr>
			<td> 881 </td>
			<td> NP6 </td>
			<td> PT6 </td>
			<td> AA </td>
		</tr>
		<tr>
			<td> 2249 </td>
			<td> NP10 </td>
			<td> PT10 </td>
			<td> AA </td>
		</tr>
		<tr>
			<td> 376 </td>
			<td> NP2 </td>
			<td> PT2 </td>
			<td> EA </td>
		</tr>
		<tr>
			<td> 647 </td>
			<td> NP5 </td>
			<td> PT5 </td>
			<td> EA </td>
		</tr>
		<tr>
			<td> 3365 </td>
			<td> NP13 </td>
			<td> PT13 </td>
			<td> EZ </td>
		</tr>
	</table>
</div>

```
ls data/*.fastq.gz
```

The FASTQ files contain the raw sequence reads for each sample. There are four lines per read: 
<ol>
	<li> Header
	<li> Sequence
	<li> Separator (usually a ‘+’)
	<li> Encoded quality value
</ol>

**Take a look at one of the FASTQ files.**

```
zless data/PT6_1.fastq.gz | head
```

Find out more about FASTQ formats at <a href="https://en.wikipedia.org/wiki/FASTQ_format/">https://en.wikipedia.org/wiki/FASTQ_format/</a>.

## 2.3 Questions

**Q1: Why is there more than one FASTQ file per sample?** *Hint: think about why there is a PT6_1.fastq.gz and a PT6_2.fastq.gz*

Now let’s move on to mapping these RNA-Seq reads to the genome using `HISAT2`.

# 3. Mapping RNA-Seq Reads to the Genome using `HISAT2`

## 3.1 Introduction

For this exercise, we have restricted the number of reads in each sample to those mapping largely to a single chromosome to reduce the mapping time. This is sufficient to illustrate the principles of differential expression analysis. 

The objectives of this part of the tutorial are to: 
<ul>
	<li> use HISAT2 to build an index from the reference genome
	<li> use HISAT2 to map RNA-Seq reads to the reference genome
</ul>

### 3.1.1 Mapping RNA-Seq reads to a genome
By this stage, you should have already performed a standard NGS quality control check on your reads to see whether there were any issues with the sample preparation or sequencing. In the interest of time, we won’t be doing that as part of this tutorial, but feel free to use the tools from earlier modules to give that a go later if you have time.

Next, we map our RNA-Seq reads to a reference genome to get context. This allows you to visually inspect your RNA-Seq data, identify contamination, novel exons and splice sites as well as giving you an overall feel for your transcriptome.

To map the RNA-Seq reads from our five samples to the reference genome, we will be using <a href="https://daehwankimlab.github.io/hisat2/">`HISAT2`</a>, a fast and sensitive splice-aware aligner. `HISAT2` compresses the genome using an indexing scheme based on the <a href="https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform">Burrows-Wheeler transform (BWT)</a> and <a href="https://en.wikipedia.org/wiki/FM-index">Ferragina-Manzini (FM)</a> index to reduce the amount of space needed to store the genome. This also makes the genome quick to search, using a whole-genome FM index to anchor each alignment and then tens of thousands local FM indexes for very rapid extensions of these alignments.

For more information, and to find the original version of _**Figure 2**_, please see the HISAT paper:
> **HISAT: a fast spliced aligner with low memory requirements Daehwan Kim.** Ben Langmead and Steven L Salzberg Nat Methods. 2015 Apr;12(4):357-60. <a href="doi:10.1038/nmeth.3317">doi:10.1038/nmeth.3317</a>

`HISAT2` is a splice-aware aligner which means it takes into account that when a read is mapped it may be split across multiple exons with (sometimes large) intronic gaps between aligned regions. As you can see in Figure 2, `HISAT2` splits read alignments into five classes based on the number of exons the read alignment is split across and the length of the anchor (longest continuously mapped portion of a split read):
<ul>
	<li> Aligns to a single exon (M) 
	<li> Alignment split across 2 exons with long anchors over 15bp (2M_gt_15) 
	<li> Alignment split across 2 exons with intermediate anchors between 8bp and 15bp (2M_8_15) 
	<li> Alignment split across 2 exons with short anchors less than 7bp (2M_1_7) 
	<li> Alignment split across more than 2 exons (gt_2M).
</ul>

`HISAT2` used the global index to place the longest continuously mapped portion of a read (anchor). This information is then used to identify the relevant local index. In most cases, `HISAT2` will only need to use a single local index to place the remaining portion of the read without having to search the rest of the genome.

For the human genome, `HISAT2` will build a single global index and 48,000 local FM indexes. Each of the local indexes represents a 64kb genomic region. The majority of human introns are significantly shorter than 64kb, so >90% of human introns fall into a single local index. Moreover, each of the local indexes overlaps its neighbour by ~1kb which means that it also has the ability to detect reads spanning multiple indexes.

There are five `HISAT2` RNA-seq read mapping categories:
<ol type="i">
	<li> M, exonic read; 
	<li> 2M_gt_15, junction reads with long, >15-bp anchors in both exons; 
	<li> 2M_8_15, junction reads with intermediate, 8- to 15-bp anchors; 
	<li> 2M_1_7, junction reads with short, 1- to 7-bp, anchors; and 
	<li> gt_2M, junction reads spanning more than two exons (Figure 2A). 
</ol>
Exoninc reads span only a single exon and represent over 60% of the read mappings in the 20 million 100-bp simulated read dataset.

## 3.2. Excercise 2

Be patient, each of the following steps will take a couple of minutes!

**Look at the usage instructions for `hisat2-build`:**

```
hisat2-build -h
```

This not only tells us the version of `HISAT2` we’re using (essential for publication methods):

>```
>HISAT2 version 2.1.0 by Daehwan Kim (infphilo@gmail.com, http://www.ccb.jhu.edu/people/infphilo
>```

But, that we also need to give hisat2-build two pieces of information:

>```
>Usage: hisat2-build [options]* <reference_in> <ht2_index_base>
>```
>These are:
><ul>
>	<li> <reference_in> location of our reference sequence file (PccAS_v3_genome.fa)
>	<li> <ht2_index_base> what we want to call our HISAT2 index files (PccAS_v3_hisat2.idx)
></ul>

**Build a `HISAT2` index for chromosome 21 of the human reference genome using `hisat2-build`:**

```
hisat2-build data/hsapien_grch38_chr21.fa outputs/hsapien_grch38_chr21_hisat2.idx
```
You can see the generated index files using:

```
ls outputs/hsapien_grch38_chr21_hisat2.idx*
```
Look at the usage for `hisat2`

```
hisat2 -h
```
Here we can see that `HISAT2` needs several parameters so that it can do the mapping:
>```
>hisat2 [options]* -x <ht2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]
>```
><ul>
>	<li> -x <ht2-idx> the prefix that we chose for our index files with hisat2-build (Pc- cAS_v3_hisat2.idx)
>	<li> {-1 <m1> -2 <m2> | -U <r>} the left (-1) and right (-2) read files for the sample (MT1_1.fastq and MT1_2.fastq respectively
>	<li> [-S <sam>] the name of the file we want to write the output alignment to (MT1.sam) as, by default, hisat2 will print the results to the terminal (stdout)

We will also be adding one more piece of information, the maximum intron length (default 500,000 bases). For this analysis, we want to set the maximum intron length to 10,000. We can do this by adding the option `--max-intronlen 10000`.

**Map the reads for the PT2 sample using `HISAT2`:**

```
 hisat2 -x outputs/hsapien_grch38_chr21_hisat2.idx -1 data/PT2_1.fastq.gz -2 data/PT2_2.fastq.gz -S outputs/PT2.sam
```

`HISAT2` has written the alignment in SAM format. This is a format which allows humans to look at our alignments. However, we need to convert the SAM file to its binary version, a BAM file. We do this for several reasons. Mainly we do it because most downstream programs require our alignments to be in BAM format and not SAM format. However, we also do it because the BAM file is smaller and so takes up less (very precious!) storage space. For more information, see the format guide: <a href="http://samtools.github.io/hts-specs/SAMv1.pdf">http://samtools.github.io/hts-specs/SAMv1.pdf</a>.

**Convert the SAM file to a BAM file:**

```
samtools view -S -o outputs/PT2.bam -b outputs/PT2.sam
```

Next we need to sort the BAM file ready for indexing. When we aligned our reads with `HISAT2`, alignments were produced in the same order as the sequences in our FASTQ files. To index the BAM file, we need the alignments ordered by their respective positions in the reference genome. We can do this using samtools sort to sort the alignments by their co-ordinates for each chromosome

**Sort the BAM file:**

```
samtools sort -o outputs/PT2_sorted.bam outputs/PT2.bam
```

Next, we need to index our BAM file. This makes searching the alignments much more efficient. It allows programs like IGV (which we will be using to visualise the alignment) to quickly get the alignments that overlap the genomic regions you’re looking at. We can do this with samtools index which will generate an index file with the extension `.bai`.

**Index the BAM file so that it can be read efficiently by `IGV`:**

```
samtools index outputs/PT2_sorted.bam
```

Now, repeat this process of mapping, converting (SAM to BAM), sorting and indexing with the reads from the NP2 sample. 

**You can run the previous steps as a single command:**

```
hisat2 -x outputs/hsapien_grch38_chr21_hisat2.idx \
	-1 data/NP2_1.fastq.gz -2 data/NP2_2.fastq.gz | \
	samtools view -S -b - | \
	samtools sort -o outputs/NP2_sorted.bam - &&
	samtools index outputs/NP2_sorted.bam
```

## 3.3. Questions
**Q1: How many index files were generated when you ran hisat2-build?** *Hint: look for the files with the .ht2 extension*

**Q2: What was the overall alignment rate for the PT2 sample to the reference genome?** *Hint: look at the the output from the HISAT2 commands*

**Q3: How does the alignment rate compare with that of the NP2 sample?** *Hint: look at the the output from the hisat2 commands*

**Q4: How many NP2 reads were not aligned to the reference genome?** *Hint: look at the the output from the hisat2 commands, you’re looking for reads (not read pairs) which have aligned 0 times (remember that one read from a pair may map even if the other doesn’t)*

The alignments generated here can be used for transcriptome quantification using tools such as feautureCounts or htseq-count and a reference transcriptome to quantify read counts at the gene or transcript level. Here we will use these alignments for visualization and then turn our attention to trancriptome alignments.

# 4. Visualising transcriptomes with `IGV`

## 4.1. Introduction
The <a href="https://software.broadinstitute.org/software/igv">Integrative Genome Viewer (IGV)</a> allows us to visualise genomic datasets.

The objectives of this part of the tutorial are:
<ul> 
	<li> load an annotation file into IGV and explore gene structure
	<li> load read alignments into IGV and inspect read alignments
</ul>

**Reference:** <a href="http://software.broadinstitute.org/software/igv/UserGuide">Online IGV User Guide</a> - see more information on all IGV features and functions.

### 4.1.1. Launch IGV and load data

```
igv &
```

This will open the IGV main window. Now, we need to tell IGV which genome we want to use. IGV has many pre-loaded genomes available so load Human (`hg38`).

### 4.1.2 Load custom gene annotation file
We not only want to see where our reads have mapped, but what genes they have mapped to. For this, we have an annotation file in GFF/GTF format. This contains a list of features, their co-ordinates and orientations which correspond to our reference genome.

<p align="center">
		<img src="images/H3ABioNet_Logo\ \(1\).png" alt="" style="width:100%">
		TEST
</p>



# 5. Transcript quantification with `Kallisto`

# 6. Differential Expression Analysis With `Sleuth`

# 7. Key Aspects of Differential Expression Analysis

# 8. Normalisation



