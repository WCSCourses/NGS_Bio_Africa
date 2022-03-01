#Module 6: Pathogen Variant Calling
##Day 1, practical session 1: Alignment to reference

### Introduction

In this module we will be analyzing the case of four Mycobacterium tuberculosis isolates (3 multi-drug resistant and 1 susceptible) from TB patients to investigate following aspects from whole genome sequence data:
1. Genetic mechanisms of resistance 
2. Identify resistance to other anti-Tb drugs 
3. Determine the genetic relatedness (pairwise SNP difference)
4. Understand their phylogenetic relationship

To answer each of these questions we will first perform the analysis on one isolate and the steps would be repeated for others as part of the assignment. 

Today we will analyze the data to answer the first two questions. Before we jump into the steps for analysis, we need to make sure we have all the required tools/softwares installed and a dataset in one place. Follow the following steps:

1. Open the Terminal and determine the current working directory
    - ```pwd```
2. Go to the folder "pathogen"
    - ```cd course_data/variant_calling/pathogen```
3. To see the contents of the folder "dataset"
    - ```ls dataset```
4. Create folder "practical" and a subfolder “MD001” within it
    - ```mkdir practical```
    - ```mkdir practical/MD001```
5. Copy the files needed for this practical session from "dataset" to "practical" folder that you created above using the following:
    - ```cp dataset/reference.fa practical/MD001/```
    - ```cp dataset/reference.gbk practical/MD001/```
6. Copying the read files:
    - ```cp dataset/MD001_R1.gz practical/MD001/```
    - ```cp dataset/MD001_R2.gz practical/MD001/```
      
     or  
      
    - ```cp dataset/MD001_*.gz practical/MD001/```  
7. To ensure the tools are installed properly the following commands when typed in the terminal must not generate any error. The command igv will open another window, please close it we will require this later in the practical.
    - ```bwa```
    - ```samtools```
    - ```picard -h```
    - ```bcftools```
    - ```snpEff -h```
    - ```igv```

### Session 1: Mapping the reads to the reference genome


1. Go to the folder "MD001" that you created above in step 4
   - ```cd practical/MD001```
2. Create index for the reference file
   - ```bwa index reference.fa```
   - The index allows the sequence aligner to more quickly find sequences in the fasta file. Like the index in a book telling you that the entry on birds starts on page 230, the genome index can tell the aligner where in the fasta file chromosome 3 is found, and how long it is. 
3. Map the reads of strain MD001 to the reference genome\
   - ```bwa mem reference.fa MD001_R1.fastq.gz MD001_R2.fastq.gz >MD001_aln.sam```
   - Note: this step might take a while to complete, so it would be worth going over the steps again to clear your doubts!!
4. Convert the alignment file " MD001_aln.sam " from SAM to BAM format
   - samtools view -O BAM -o MD001_aln.bam MD001_aln.sam
   - You can check the ".sam” and “.bam” file you just created (hint: use ls -lh) to observe how much of a difference it makes by converting .sam to .bam files
5. Sort the alignment file "MD001_aln.bam"
   - ```samtools sort -T temp -O bam -o MD001_aln_sorted.bam MD001_aln.bam```
6. Index the sorted bam file
   - ```samtools index MD001_aln_sorted.bam```
7. Mark the duplicates (these are potential optical duplicates, please refer to the QC material for better understanding)
   - ```picard MarkDuplicates I=MD001_aln_sorted.bam O=MD001_aln_markdup.bam M=metrics.txt```
   - In order to check how many reads were marked as duplicates we can use the following command:
   - ```grep -A 2 “^##METRICS” metrics.txt ```
   - Here, `-A` option allows us to print number of lines (2 in above case) after the line with the matching string `(“##METRICS)`

**Q1.1:** How many read pairs were assigned as duplicates?

**Q1.2:** What proportion of the mapped sequence was marked as duplicates?

8. Index the file created above "MD001_aln_markdup.bam"
   - ```samtools index MD001_aln_markdup.bam```
9. Generate the statistics of mapping using the command
   - ```samtools stats MD001_aln_markdup.bam >MD001_bamstats.txt```
   - ```grep “^SN” MD001_bamstats.txt > MD001stats.txt```

**Q1.3:** What is the total number of mapped reads?

**Q1.4:** What is the total number of unmapped reads?

**Q1.5:** What is the total number of mapped and properly paired reads?

**Q1.6:** What is the average insert size?

**Q1.7:** What is the percentage of reads properly paired?

## Session 2: Variant calling

1. Extract all the variants from the alignment file "MD001_aln_markdup.bam"
   - ```samtools mpileup -t DP -Bug -m 4 -f reference.fa MD001_aln_markdup.bam >MD001_variants.bcf```
   - Note: You might see a few warnings, you can ignore them
2. Convert bcf to vcf file
   - ```bcftools call -mv -O v -o MD001_variants.vcf MD001_variants.bcf```
3. To see the first 100 lines of the "MD001_variants.vcf" file
   - ```head -n100 MD001_variants.vcf```

**Q2.1:** At what position is the first variant in the unfiltered vcf file for MD001?

4. Filtering the variants
   - At this stage the variants include both SNPs and short Indels identified from the alignment of reads. The following process describes how we can identify high quality SNPs (only) from the file applying different metrics
5. Filter only SNPs from the "MD001_variants.vcf" file
   - ```bcftools filter -i 'type="snp"' -g10 -G10 MD001_variants.vcf -o MD001_SNPs.vcf```
6. Filter the SNPs with base quality (QUAL) >=50, MQ >=30 and read depth (DP) >5
   - ```bcftools filter -i 'type="snp" && QUAL>=50 && FORMAT/DP>5 && MQ>=30' -g10 -G10 MD001_variants.vcf -o MD001_SNPs_filtered_try1.vcf```
   - Note: There is no space between the two “&&” symbols.
7. Filter all homozygous SNPs (alternate base ratio >= 0.80)
   - ```bcftools filter -i 'type="snp" && QUAL>=50 && FORMAT/DP>5 && MQ>=30 && DP4[2]/(DP4[2]+DP4[0])>=0.80 && DP4[3]/(DP4[3]+DP4[1])>=0.80' -g10 -G10 MD001_variants.vcf -o MD001_SNPs_filtered.vcf```

**Q2.2:** What does the DP4 value represent?

8. Annotating variants
   - Annotating the variants allows us to investigate what effect of the mutations are, and which have mutations in drug resistance related genes.
9. Find the correct database for snpEff
   - ```snpEff databases | grep h37rv```
10. 2.10 Rename the chromosome
    - ```sed ‘s|NC_000962.3|Chromosome|’ MD001_SNPs_filtered.vcf > MD001_SNPs_filtered_renamed.vcf```
11. Run the annotation 
    - ```snpEff \
        Mycobacterium_tuberculosis_h37rv \
        -csvStats MD001_snpEff.csv \
        -v \
        MD001_SNPs_filtered_renamed.vcf \
        > MD001_SNPs_filtered_snpEff.ann.vcf
      ```
12. Identify for variants in candidate drug resistance genes 
    - ```cat MD001_SNPs_filtered_snpEff.ann.vcf | grep MODERATE\|rpoB```

**Q2.5:** How many HIGH effect variants were there for sample MD001?

**Q2.6:** What was the TS TV ratio? (Look in the MD001_snpEff.csv summary file)

**Q2.7:** Are there any other mutations in resistance related genes? 






