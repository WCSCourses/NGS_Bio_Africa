## 1. RNA-Seq Expression Analysis

### 1.1. Introduction

RNA sequencing (RNA-Seq) is a high-throughput method used to profile the transcriptome, quan- tify gene expression and discover novel RNA molecules or splice variants. This tutorial uses RNA sequencing of human cancer samples to walk you through transcriptome alignment, visualisation, simple quality control checks and shows you how to profile transcriptomic differences by identifying differentially expressed genes.

For an introduction to RNA-Seq principles and best practices see:

> **A survey of best practices for RNA-Seq data analysis**  
Ana Conesa, Pedro Madrigal, Sonia Tarazona, David Gomez-Cabrero, Alejandra Cervera, Andrew McPherson, Michał Wojciech Szcześniak, Daniel J. Gaffney, Laura L. Elo, Xue- gong Zhang and Ali Mortazavi Genome Biol. 2016 Jan 26;17:13 doi:10.1186/s13059-016-0881-8

### 1.2. Learning Outcomes

By the end of this tutorial you can expect to be able to:
<ul>
	<li> Align RNA-Seq reads to a reference genome and a transcriptome
	<li> Visualise transcription data using standard tools
	<li> Perform QC of NGS transcriptomic data
	<li> Quantify the expression values of your transcripts using standard tools
</ul>

### 1.3. Practical Outline

This tutorial comprises the following sections:  
<ol>
	<li>Introducing the tutorial dataset  
	<li>Mapping RNA-Seq reads to the genome with HISAT2  
	<li>Visualising transcriptomes with IGV  
	<li>Transcript quantification with Kallisto  
	<li>Identifying differentially expressed genes with Sleuth  
	<li>Key aspects of differential expression analysis
</ol>

### 1.4. Authors

This tutorial was developed by Victoria Offord and Adam Reid and adapted for use here by Nyasha Chambwe.

### 1.5. Prerequisites

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

## 2. Introducing The Tutorial Dataset

## 3. Mapping RNA-Seq Reads to the Genome using `HISAT2`

## 4. Visualising transcriptomes with IGV

## 5. Transcript quantification with Kallisto

## 6. Differential Expression Analysis With Sleuth

## 7. Key Aspects of Differential Expression Analysis

## 8. Normalisation



