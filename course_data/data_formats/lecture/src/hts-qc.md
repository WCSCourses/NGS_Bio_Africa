% only for printing
%   STYLE:\setbeamertemplate{footline}{}
% followed by
%   subset-slides

\includegraphics[width=1.0\linewidth]{img/wtsi.png}
\vskip2em
\centerline{\large HTS data formats and Quality Control}\smallskip
\centerline{petr.danecek@sanger.ac.uk}
\vskip10em
\begin{textblock}{4}(0,7)
\includegraphics[width=2.5\TPHorizModule]{img/wtsi-logo.png}
\end{textblock}




# Data Formats
FASTQ

* Unaligned read sequences with base qualities

SAM/BAM

* Unaligned or aligned reads
* Text and binary formats

CRAM

* Better compression than BAM

VCF/BCF

* Flexible variant call format
* Arbitrary types of sequence variation
* SNPs, indels, structural variations

\vskip2em

Specifications maintained by the Global Alliance for Genomics and Health

\begin{textblock}{2}(7.5,1.2)
\includegraphics[width=1.4\TPHorizModule]{img/fastq-sam-vcf.pdf}
\end{textblock}




# FASTQ {.fragile}
* Simple format for raw unaligned sequencing reads
* Extension to the FASTA file format
* Sequence and an associated per base quality score

\begin{Verbatim}[fontsize=\footnotesize,xleftmargin=3em,commandchars=\\\{\}]
\textbf{\color{darkgold}@ERR007731.739 IL16_2979:6:1:9:1684/1}
\textbf{\color{darkgold}CTTGACGACTTGAAAAATGACGAAATCACTAAAAAACGTGAAAAATGAGAAATG}
\textbf{\color{darkgold}+}
\textbf{\color{darkgold}BBCBCBBBBBBBABBABBBBBBBABBBBBBBBBBBBBBABAAAABBBBB=@>B}
@ERR007731.740 IL16_2979:6:1:9:1419/1
AAAAAAAAAGATGTCATCAGCACATCAGAAAAGAAGGCAACTTTAAAACTTTTC
+
BBABBBABABAABABABBABBBAAA>@B@BBAA@4AAA>.>BAA@779:AAA@A
\end{Verbatim}

* Quality encoded in ASCII characters with decimal codes 33-126

    * ASCII code of "A" is 65, the corresponding quality is Q$= 65 - 33 = 32$

    \vskip0em

    * Phred quality score: $P = 10^{-Q/10}$
    \hbox{\small\hskip4.5em\tt
        perl -e 'printf "\%d\textbackslash n",ord("A")-33;'
    }

* Beware: multiple quality scores were in use!
    * Sanger, Solexa, Illumina 1.3+

* Paired-end sequencing produces two FASTQ files 


\pause
\vonly[2-]{{\color{red}Q: What is the probability of a sequencing error if the quality is "?" (ASCII code 63)}}

\pause
{\color{darkgreen}A: Q=30, one error in 1000 bases}

\note[item]{Base calling quality: originally determined from Sanger sequencing
    metrics such as peak resolution and shape.
    Today each sequencing chemistry has relevant parameters, these are
    are compared to a large empirical data set of known accuracy.
    The parameters are: intensity profiles, signal-to-noise ratio, possibly others.
}



# SAM / BAM
SAM (Sequence Alignment/Map) format

* Unified format for storing read alignments to a reference genome
* Developed by the 1000 Genomes Project group (2009)

\pause

* One record (a single DNA fragment alignment) per line describing alignment between fragment and reference
* 11 fixed columns + optional key:type:value tuples

\vskip0.5em
\centerline{\includegraphics[width=7\TPHorizModule]{img/sam-record.pdf}}
\vskip0.5em

\pause 

\vskip1em

Note that BAM can contain

* unmapped reads
* multiple alignments of the same read
* supplementary (chimeric) reads


# SAM {.fragile}

SAM fields
\begingroup
\footnotesize
\begin{tabular}{l l l}
1      &     QNAME     &      Query NAME of the read or the read pair               \\
2      &     FLAG      &      Bitwise FLAG (pairing, strand, mate strand, etc.)     \\
3      &     RNAME     &      Reference sequence NAME                               \\
4      &     POS       &      1-Based leftmost POSition of clipped alignment        \\
5      &     MAPQ      &      MAPping Quality (Phred-scaled)                        \\
6      &     CIGAR     &      Extended CIGAR string (operations: MIDNSHPX=)         \\
7      &     MRNM      &      Mate Reference NaMe ('=' if same as RNAME)            \\
8      &     MPOS      &      1-Based leftmost Mate POSition                        \\
9      &     ISIZE     &      Inferred Insert SIZE                                  \\
10     &     SEQ       &      Query SEQuence on the same strand as the reference    \\
11     &     QUAL      &      Query QUALity (ASCII-33=Phred base quality)           \\
12-    &     OTHER     &      Optional fields 
\end{tabular}
\vskip1em
\endgroup

\hbox{\footnotesize\vbox{\begin{verbatim}
$ samtools view -h file.bam | less

@HD VN:1.0  GO:none SO:coordinate
@SQ SN:1    LN:249250621    UR:hs37d5.fa.gz AS:NCBI37   M5:1b22b98cdeb4a9304cb5d48026a85128 SP:Human
@SQ SN:2    LN:243199373    UR:hs37d5.fa.gz AS:NCBI37   M5:a0d9851da00400dec1098a9255ac712e SP:Human
@RG ID:1    PL:ILLUMINA PU:13350_1  LB:13350_1  SM:13350_1  CN:SC
@PG ID:bwa  PN:bwa  VN:0.7.10-r806    CL:bwa mem hs37d5.fa.gz 13350_1_1.fq 13350_1_1.fq
1:2203:10256:56986  97  1   9998    0   106M45S =  10335    0 \
    CCATAACCCTAACCCTAACCCTAACCATAGCCCTAACCCTAACCCTAACCCTAACCCT[...]CAAACCCACCCCCAAACCCAAAACCTCACCAC \
    FFFFFJJJJJJJJFJJJJFJAJJJJJ-JJAAAJFJJFFJJF<FJJFFJJJJFJJJJFF[...]<---F-----A7-J-<J-A--77AF---J7-- \
    MD:Z:1G24C2A76  PG:Z:MarkDuplicates RG:Z:1  NM:i:3  MQ:i:0  AS:i:94 XS:i:94
\end{verbatim}}}



# CIGAR string {.fragile}

Compact representation of sequence alignment

\begingroup
\footnotesize
\begin{tabular}{l l l}
M   &   alignment match or mismatch                                 \\
=   &   sequence match                                              \\
X   &   sequence mismatch                                           \\
I   &   insertion to the reference                                  \\
D   &   deletion from the reference                                 \\
S   &   soft clipping (clipped sequences present in SEQ)            \\
H   &   hard clipping (clipped sequences NOT present in SEQ)        \\
N   &   skipped region from the reference                           \\
P   &   padding (silent deletion from padded reference) 
\end{tabular}
\vskip1em
\endgroup

Examples:
\begin{Verbatim}[commandchars=\\\{\},fontsize=\footnotesize,xleftmargin=2em]
Ref:   ACGTACGTACGTACGT
Read:  ACGT----ACGTACGA
Cigar: 4M 4D 8M

Ref:   ACGT----ACGTACGT
Read:  ACGTACGTACGTACGT
Cigar: 4M 4I 8M

Ref:   ACTCAGTG--GT
Read:  ACGCA-TGCAGTtagacgt
Cigar: 5M 1D 2M 2I 2M 7S

\vonly[2-]{Ref:   tgtcgtcACGCATG---CAGT}
\vonly[2-]{Read:         AAGCATGCGGCAGTacgtacg}
\vonly[2-]{\color{red}Cigar:}\vonly[3]{  7H 7M 3I 4M 7S}
\end{Verbatim}



# Flags

\vskip0.5em 

\footnotesize
\begin{tabular}{l l l l}
Hex   & Dec  &  Flag           & Description \\
\hline \\
0x1   & 1    &  PAIRED         & paired-end (or multiple-segment) sequencing technology    \\ 
0x2   & 2    &  PROPER\_PAIR   & each segment properly aligned according to the aligner    \\ 
0x4   & 4    &  UNMAP          & segment unmapped                                          \\ 
0x8   & 8    &  MUNMAP         & next segment in the template unmapped                     \\ 
0x10  & 16   &  REVERSE        & SEQ is reverse complemented                               \\ 
0x20  & 32   &  MREVERSE       & SEQ of the next segment in the template is reversed       \\ 
0x40  & 64   &  READ1          & the first segment in the template                         \\ 
0x80  & 128  &  READ2          & the last segment in the template                          \\ 
0x100 & 256  &  SECONDARY      & secondary alignment                                       \\ 
0x200 & 512  &  QCFAIL         & not passing quality controls                              \\ 
0x400 & 1024 &  DUP            & PCR or optical duplicate                                  \\ 
0x800 & 2048 &  SUPPLEMENTARY  & supplementary alignment                                   \\ 
\end{tabular}
\normalsize

\vskip0.5em 

Bit operations made easy

* python
\hbox{\hskip1em\vbox{
\footnotesize
    0x1 | 0x2 | 0x20 | 0x80     ..  163 \\
    bin(163)   ..  10100011
\normalsize}}

* samtools flags
\hbox{\hskip1em\vbox{
\footnotesize
    0xa3    163 PAIRED,PROPER\_PAIR,MREVERSE,READ2
\normalsize}}

\vskip5em

\begin{textblock}{8}(0,7.2)
\vonly[2-4]{{\color{red}Q: What is the flag of the first read if both reads are mapped?}} \vskip0em
\vonly[3]{{\color{darkgreen}\small A: PAIRED,PROPER\_PAIR,READ1 .. \normalsize}}
\vonly[4]{{\color{darkgreen}\small A: PAIRED,PROPER\_PAIR,READ1 .. 1+2+64=67\normalsize}}
\end{textblock}
\begin{textblock}{8}(0,7.2)
\vonly[5-6]{{\color{red}Q: What is the meaning of the flag 65?}} \vskip0em
\vonly[6]{{\color{darkgreen}\small A: 65=1+64 .. PAIRED,READ1\normalsize}}
\end{textblock}




# Optional tags {.fragile}

Each lane has a unique RG tag that contains meta-data for the lane

RG tags

* ID: SRR/ERR number

* PL: Sequencing platform

* PU: Run name

* LB: Library name

* PI: Insert fragment size

* SM: Individual

* CN: Sequencing center

\vskip10em

\begin{textblock}{5}(5,1.5)
\includegraphics[width=4\TPHorizModule]{img/pipeline.pdf}
\end{textblock}



# BAM

BAM (Binary Alignment/Map) format

* Binary version of SAM
* Developed for fast processing and random access
    * BGZF (Block GZIP) compression for indexing

Key features

* Can store alignments from most mappers
* Supports multiple sequencing technologies
* Supports indexing for quick retrieval/viewing
* Compact size (e.g. 112Gbp Illumina = 116GB disk space)
* Reads can be grouped into logical groups e.g. lanes, libraries, samples
* Widely supported by variant calling packages and viewers




# Reference based Compression

BAM files are too large

* ~1.5-2 bytes per base pair

Increases in disk capacity are being far outstripped by sequencing technologies

\begin{textblock}{9}(0,2.5)
\vonly[1]{\centerline{\includegraphics[width=8\TPHorizModule]{img/omics-data-size.png}}}
\end{textblock}
\begin{textblock}{9}(0,8.2)
\footnotesize 
\vonly[1]{\centerline{Zachary D. Stephens, \textit{et al}, Big Data: Astronomical or Genomical? DOI: 10.1371/journal.pbio.1002195}}
\end{textblock}

\pause

BAM stores all of the data

* Every read base
* Every base quality
* Using a single conventional compression technique for all types of data

\vskip14em

\begin{textblock}{9}(0,4.5)
\vonly[2]{\centerline{\includegraphics[width=7\TPHorizModule]{img/ref-based-compression1.pdf}}}
\vonly[3]{\centerline{\includegraphics[width=7\TPHorizModule]{img/ref-based-compression2.pdf}}}
\end{textblock}



# CRAM

Three important concepts

* Reference based compression
* Controlled loss of quality information
* Different compression methods to suit the type of data, e.g. base qualities vs. metadata vs. extra tags

In lossless mode: 60% of BAM size

Archives and sequencing centers moving from BAM to CRAM

* Support for CRAM added to Samtools/HTSlib in 2014
* Soon to be available in Picard/GATK

\vskip1em
\centerline{\includegraphics[width=6\TPHorizModule]{img/cram.pdf}}



# VCF: Variant Call Format

File format for storing variation data

* Tab-delimited text, parsable by standard UNIX commands

* Flexible and user-extensible

* Compressed with BGZF (bgzip), indexed with TBI or CSI (tabix)

\vskip1em
\centerline{\includegraphics[width=10\TPHorizModule]{img/vcf.png}}
\vskip3em


# VCF / BCF

VCFs can be very big

* compressed VCF with 3781 samples, human data: 
    * 54 GB for chromosome 1
    * 680 GB whole genome

VCFs can be slow to parse

* text conversion is slow

* main bottleneck: FORMAT fields

\vskip0.5em
\centerline{\includegraphics[width=10\TPHorizModule]{img/vcf-to-bcf.png}}
\vskip1em

BCF 

* binary representation of VCF

* fields rearranged for fast access

\vskip0.5em
\centerline{\includegraphics[width=10\TPHorizModule]{img/vcf-to-bcf2.png}}


# gVCF

Often it is not enough not know *variant* sites only

* was a site dropped because of a reference call or because of missing data?
* we need evidence for both variant and non-variant positions in the genome

\vskip0.5em

gVCF

* blocks of reference-only sites can be represented in a single record using the INFO/END tag

* symbolic alleles <*> for incremental calling
    * raw, "callable" gVCF
    * calculate genotype likelihoods only once (an expensive step)
    * then call incrementally as more samples come in

\vskip0.5em

\centerline{\includegraphics[width=1\textwidth]{img/gvcf.pdf}}




# Optimizing variant calls for speed

\vskip2em

\centerline{\includegraphics[width=1\textwidth]{img/tomahawk.pdf}}

\vskip3em

New TWK format by Marcus Klarqvist (under development)

- BCF still too slow for querying hundreds of thousands and millions of samples 

- bigger but 100x faster for certain operations on GTs



# Custom formats for custom tasks

\centerline{
\vbox{\hsize=0.5\textwidth
    \centerline{Region in LD}
    \vskip0.5em
    \centerline{\includegraphics[width=0.5\textwidth]{img/region-in-LD.png}}
}
\hskip1em
\vbox{\hsize=0.5\textwidth
    \centerline{Random region}
    \vskip0.5em
    \centerline{\includegraphics[width=0.5\textwidth]{img/random-region.png}}
}}
\vskip1em
\centerline{
    \includegraphics[width=0.4\textwidth]{img/ld-scaffold.png}
    \hskip2em
    \includegraphics[width=0.4\textwidth]{img/ld-scaffold-graph.png}
}

\begin{textblock}{9}(0,8.5)
\footnotesize 
\centerline{Method and figures by Marcus Klarqvist}
\end{textblock}





# Global Alliance for Genomics and Health

International coalition dedicated to improving human health 

Mission

* establish a common framework to enable sharing of genomic and clinical data

Working groups

* clinical
* regulatory and ethics
* security
* data

Data working group

* beacon project \footnotesize .. test the willingness of international sites to share genetic data \normalsize
* BRCA challenge \footnotesize .. advance understanding of the genetic basis of breast and other cancers \normalsize
* matchmaker exchange \footnotesize .. locate data on rare phenotypes or genotypes \normalsize
* reference variation \footnotesize .. describe how genomes differ so researchers can assemble and interpret them \normalsize
* benchmarking \footnotesize .. develop variant calling benchmark toolkits for germline, cancer, and transcripts \normalsize
* file formats \footnotesize .. CRAM, SAM/BAM, VCF/BCF \normalsize

File formats

* http://samtools.github.io/hts-specs/

\begin{textblock}{4}(4.5,2.5)
\includegraphics[width=4\TPHorizModule]{img/globalalliance.jpg}
\end{textblock}


%   # Coffee break and questions
%   
%   \vskip7em
%   \centerline{\includegraphics[width=2\TPHorizModule]{img/coffee-break.pdf}}



# Quality Control

Biases in sequencing

* Base calling accuracy

* Read cycle vs. base content

* GC vs. depth

* Indel ratio

\vskip0.5em

Biases in mapping

\vskip1em

Genotype checking

* Sample swaps

* Contaminations



# Base quality

Sequencing by synthesis: dephasing

* growing sequences in a cluster gradually desynchronize
* error rate increases with read length

\vskip1em

Calculate the average quality at each position across all reads

\vskip0.5em
\centerline{\includegraphics[width=5\TPHorizModule]{img/read-qual.png}}

\vskip0.5em
\centerline{\hbox{
\footnotesize
\begin{tabular}{l l l}
Quality   &   Probability of error  &  Call Accuracy   \\ 
10 (Q10)  &        1 in 10          &  90\%            \\ 
20 (Q20)  &        1 in 100         &  99\%            \\ 
30 (Q30)  &        1 in 1000        &  99.9\%           \\
40 (Q40)  &        1 in 10000       &  99.99\%                
\end{tabular}
\normalsize
}}


# Base calling errors

\centerline{\includegraphics[width=9\TPHorizModule]{img/base-calling-biases.jpg}}

\begin{textblock}{9}(0,8.5)
\footnotesize 
\centerline{Base-calling for next-generation sequencing platforms, doi: 10.1093/bib/bbq077}
\end{textblock}

<!-- signal decay caused by the DNA loss due to primer-template melting, digestion by enzymatic impurities, DNA dissociation and misincorporation -->



# Base quality

\centerline{
\includegraphics[height=5\TPHorizModule]{{img/read-qual.pass}.png}\hskip5em
\includegraphics[height=5\TPHorizModule]{{img/read-qual.fail}.png}
}
\vonly[2]{
\begin{textblock}{4}(2,7.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,7.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}



# Mismatches per cycle

Mismatches in aligned reads (requires reference sequence)

* detect cycle-specific errors
* base qualities are informative!


\vskip1em
\centerline{
\includegraphics[width=5\TPHorizModule]{{img/mismatch-per-cycle.pass}.png}
\includegraphics[width=5\TPHorizModule]{{img/mismatch-per-cycle.fail}.png}
}

\vonly[2]{
\begin{textblock}{4}(2,7.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,7.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}



# GC bias

GC- and AT-rich regions are more difficult to amplify

* compare the GC content against the expected distribution (reference sequence)

\vskip2em
\centerline{
\includegraphics[width=10\TPHorizModule]{img/gc-content.png}
}
\vskip3em

\vonly[2-3]{
\begin{textblock}{4}(2,7)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,7)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}
\vonly[3]{
\begin{textblock}{4}(4.8,2.1)
\includegraphics[width=4.5\TPHorizModule]{{img/gc-content.fail}.png}
\end{textblock}
}



# GC content by cycle

Was the adapter sequence trimmed?

\vskip1em
\centerline{
\includegraphics[width=4.9\TPHorizModule]{{img/acgt-per-cycle.pass}.png}
\includegraphics[width=5.1\TPHorizModule]{{img/acgt-per-cycle.fail}.png}
}
\vskip3em

\vonly[2-]{
\begin{textblock}{4}(2,6.5)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,6.5)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}


# Fragment size

Paired-end sequencing: the size of DNA fragments matters

\vskip1em
\vonly[1-2]{
\centerline{
\includegraphics[width=5\TPHorizModule]{{img/insert-size.pass}.png}
\includegraphics[width=5\TPHorizModule]{{img/insert-size.fail}.png}
}}
\vskip3em

\vonly[2]{
\begin{textblock}{4}(2,6.5)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,6.5)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}



# Quiz


\vskip1em
\centerline{\includegraphics[width=5\TPHorizModule]{{img/insert-size.fail2}.png}\hskip2em}
\vskip1em
\centerline{This is 100bp paired-end sequencing. Can you spot any problems??}
\vskip3em
\pause
\centerline{\includegraphics[width=6\TPHorizModule]{img/fragment-size-quizz.pdf}}



# Insertions / Deletions per cycle

False indels

* air bubbles in the flow cell can manifest as false indels

\vskip1em
\centerline{
\includegraphics[width=5\TPHorizModule]{{img/indels-per-cycle.pass}.png}
\includegraphics[width=5\TPHorizModule]{{img/indels-per-cycle.fail}.png}
}
\vskip3em

\vonly[2]{
\begin{textblock}{4}(2,7)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,7)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}



# Auto QC tests

A suggestion for human data:
\vskip1em
\small
\begin{tabular}{l r}
Minimum number of mapped bases & 90\% \\
Maximum error rate & 0.02\% \\
Maximum number of duplicate reads & 5\% \\
Minimum number of mapped reads which are properly paired & 80\% \\
Maximum number of duplicated bases due to overlapping read pairs & 4\% \\
Maximum in/del ratio & 0.82 \\
Minimum in/del ratio & 0.68 \\
Maximimum indels per cycle, factor above median & 8 \\
Minimum number of reads within 25\% of the main peak & 80\% \\
\end{tabular}

\vskip2em
\centerline{\includegraphics[width=4\TPHorizModule]{img/isize.pdf}}     <!-- isize.py -->


# Detecting sample swaps

Check the identity against a known set of variants

\vskip1em
\centerline{
\includegraphics[width=4.5\TPHorizModule]{{img/gtcheck.dense}.pdf}\hskip1em
\includegraphics[width=4.5\TPHorizModule]{{img/gtcheck.fail}.pdf}
}
\vskip3em

\vonly[2]{
\begin{textblock}{4}(2,6.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-yes.pdf}
\end{textblock}
\begin{textblock}{4}(7,6.3)
\includegraphics[width=0.5\TPHorizModule]{img/tickmark-no.pdf}
\end{textblock}
}



# How many markers are necessary?

\vskip1em
Number of sites required to identify non-related human samples

\vskip1em
\centerline{
\includegraphics[width=4.5\TPHorizModule]{{img/gtcheck.dense}.pdf}\hskip1em
\includegraphics[width=4.5\TPHorizModule]{{img/gtcheck.sparse}.pdf}
}
\centerline{
\includegraphics[width=4.5\TPHorizModule]{{img/gtcheck.exsparse}.pdf}
}


# Software

Software used to produce graphs in these slides

* `samtools stats` and `plot-bamstats`

* `bcftools gtcheck` 

* `matplotlib`


%   # xxx
%   
%   Pipelining
%   
%   http://seqanswers.com
%   http://www.cbs.dtu.dk/courses/27626/Exercises/BAM-postprocessing.php
%


% Schwartz: Detection and Removal of Biases in the Analysis of Next- Generation Sequencing Reads 
%
%   # Biases in Next Generation Sequencing data
%   
%   * nucleotide per cycle bias
%       * mostly in RNA-seq, sometimes in ChIP-seq
%       * cannot be attributed to biased PCR-amplification
%       * partial explanation: random hexamer priming during reverse transcription?
%       - more references in the paper above: 14,15,16,17
%   
%   dephasing
%   
%   The illumina platform uses a so called sequencing by synthesis process. Bases are added one at a time and the consensus is determined in a cluster of identical sequences.
%   
%   The source of errors can be numerous, here is one review that discusses the issues in more detail:
%   
%   The challenges of sequencing by synthesis, Nature Biotech, 2009
%   
%   In a nuthsell a short answer to the best of my understanding is this: Not
%   all sequences in a cluster will grow at the same rate, this will slowly
%   lead to a desynchronization as the errors accumulate. This is why the
%   quality dips towards the end.



