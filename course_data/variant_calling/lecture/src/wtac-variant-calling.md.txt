% only for printing
%   STYLE:\setbeamertemplate{footline}{}
% followed by
%   subset-slides


% https://docs.google.com/presentation/d/1Log3pInsXvh-vKJQsQiAWRvq6GoMuZtEmAiUAo5tem8/edit?ts=57ce8b1d#slide=id.p
% https://docs.google.com/document/d/1J3RWKQKU2x3th1QI4e14lK5ATj6Seee-4QyxMkxQ-_g/edit?ts=57ce8b2b

% synthetic diploid benchmark
% !! https://www.biorxiv.org/content/early/2017/11/22/223297


\includegraphics[width=1.0\linewidth]{img/wtsi.png}
\vskip2em
\centerline{\large Variant Calling - SNPs and short indels}\smallskip
\centerline{petr.danecek@sanger.ac.uk}
\vskip10em
\begin{textblock}{4}(0,7)
\includegraphics[width=2.5\TPHorizModule]{img/wtsi-logo.png}
\end{textblock}




# HTS workflow

Library preparation

* DNA extraction
* fragmentation
* adapter ligation
* amplification

Sequencing

* base calling
* de-multiplexing

Data processing

* read mapping
* variant calling
* variant filtering

Analysis

* Variant annotation
* $\ldots$

\begin{textblock}{2}(7,0.8)
\includegraphics[width=1.5\TPHorizModule]{img/hts-workflow.pdf}
\end{textblock}


# Variant types

SNPs/SNVs \hskip0.5em \small$\ldots$ Single Nucleotide Polymorphism/Variation \normalsize

\vskip0.8em
\hskip4em\includegraphics[height=0.5\TPHorizModule]{img/what-is-snp.pdf}
\vskip1em

MNPs \hskip0.5em \small$\ldots$ Multi-Nucleotide Polymorphism \normalsize

\vskip0.8em
\hskip4em\includegraphics[height=0.5\TPHorizModule]{img/what-is-mnp.pdf}
\vskip1em

Indels \hskip0.5em \small$\ldots$ short insertions and deletions \normalsize

\vskip0.8em
\hskip4em\includegraphics[height=0.5\TPHorizModule]{img/what-is-indel.pdf}
\vskip1em

SVs \hskip0.5em \small$\ldots$ Structural Variation \normalsize

\vskip0.8em
\hskip4em\includegraphics[height=1.8\TPHorizModule]{img/what-is-sv.pdf}




# Some terminology

The goal is to determine the genotype at each position in the genome 

\vskip1em

Genotype

* in the broad sense ... genetic makeup of an organism
* in the narrow sense ... the combination of alleles at a position

\vskip1em

Reference and alternate alleles - R and A

\vskip1em

Diploid organism

* two chromosomal copies, three possible genotypes
    * RR .. homozygous reference genotype
    * RA .. heterozygous
    * AA .. homozygous alternate

\vskip1em
\centerline{\includegraphics[width=6.5\TPHorizModule]{img/diploid-example.pdf}}


# Germline vs somatic mutation

Germline mutation

* heritable variation in the germ cells

\vskip0.5em

Somatic mutation

* variation in non-germline tissue, tumors...

\pause
\vskip1em

Germline variant calling 

* expect the following fractions of alternate alleles in the pileup:\small
\vskip0.1em\hskip2em 0.0 for RR genotype (plus sequencing errors)
\vskip0.1em\hskip2em 1.0 for AA (plus sequencing errors)
\vskip0.1em\hskip2em 0.5 for RA (random variation of binomial sampling)

\normalsize

Somatic

* any fraction of alt AF possible - subclonal variation, admixture of normal cells in tumor sample 

\vskip1em
\centerline{\includegraphics[width=5.5\TPHorizModule]{img/pileup.pdf}\hskip2em}


# Naive variant calling

Use fixed allele frequency threshold to determine the genotype

\vskip1em
\vonly[1]{\centerline{\includegraphics[width=6\TPHorizModule]{img/naive-calling-0.pdf}\hskip2em}}
\vonly[2]{\centerline{\includegraphics[width=6\TPHorizModule]{img/naive-calling-1.pdf}\hskip2em}}
\vonly[3-]{\centerline{\includegraphics[width=6\TPHorizModule]{img/naive-calling-2.pdf}\hskip2em}}

\vskip2em

\pause

\noindent 1) Filter base calls by quality

\hskip2em e.g. ignore bases Q<20

\vskip1em
\pause

\noindent 2) Filter reads with low mapping quality

\pause
\vskip1em

Problems:

* undercalls hets in low-coverage data
* throws away information due to hard quality thresholds
* gives no measure of confidence


\vskip8em

\vonly[2]{
\begin{textblock}{1}(4.2,3.8)
\includegraphics[width=5\TPHorizModule]{img/read-qual.png}
\end{textblock}
}
\vonly[2]{
\begin{textblock}{4}(-0.1,5.45)
\small
\hskip0.7em {\bf Phred quality score} \\
\hskip0.7em $\mbox{Q} = -10 \log_{10} P_{\mbox{err}}$
\vskip1em
\footnotesize
\begin{tabular}{l l l}
Quality   &   Error probability &  Accuracy   \\ 
10 (Q10)  &        1 in 10          &  90\%            \\ 
20 (Q20)  &        1 in 100         &  99\%            \\ 
30 (Q30)  &        1 in 1000        &  99.9\%           \\
40 (Q40)  &        1 in 10000       &  99.99\%                
\end{tabular}
\end{textblock}}

\vonly[1,3,4]{
\begin{textblock}{4}(5,4.2)
\centerline{\hbox{
\small
\begin{tabular}{|l l|}
\hline 
\rowcolor[gray]{0.8}
\hskip0.1em alt AF        &   \hskip2em genotype \vbox to1.2em{}\\
\hline 
$[0,0.2)$     &   RR .. homozygous reference \vbox to1.2em{}\\ 
$[0.2,0.8]$   &   RA .. herezogyous \\ 
$(0.8,1]$     &   AA .. homozygous variant  \\
\hline 
\end{tabular}
\normalsize
}}
\end{textblock}}




# Real life calling models

More sophisticated models apply a statistical framework
\vskip1em

\centerline{\includegraphics[width=3.5\TPHorizModule]{img/bayes-theorem-pdg.pdf}}

to determine:
\vskip1em

1. the most likely genotype $g\in\{\mbox{RR},\mbox{RA},\mbox{AA}\}$ given the observed data $D$
$$g = \operatornamewithlimits{argmax}_G P(G|D) \hskip12em$$

\vskip1em

2. and the genotype quality
$$ Q = -10 \log_{10} [1 - P(G|D)] \hskip12em$$

\begin{textblock}{1}(6,5.7)
\hskip2em\includegraphics[width=3\TPHorizModule]{img/phred-score-pdg.pdf}
\end{textblock}



# Important terms you may encounter

\definecolor{mblue}{HTML}{337ab7}{\color{mblue}Genotype likelihoods}

* which of the three genotypes RR, RA, AA is the data most consistent with?

* calculated from the alignments, the basis for calling

* takes into account:
    * base calling errors
    * mapping errors
    * statistical fluctuations of random sampling
    * local indel realignment (base alignment quality, BAQ)

\vskip1em

\definecolor{morange}{HTML}{F0AD4E}{\color{morange}Prior probability}

* how likely it is to encounter a variant base in the genome?

* some assumptions are made
    * allele frequencies are in Hardy-Weinberg equilibrium
    \vskip0.3em\hskip1em \footnotesize $P(\mbox{RA}) = 2 f (1-f)$, $P(\mbox{RR}) = (1-f)^2$, $P(\mbox{AA}) = f^2$

* can take into account genetic diversity in a population

\vskip1em

\centerline{\includegraphics[width=3\TPHorizModule]{img/bayes-theorem-pdg2.pdf}}



# Variant calling example

Inputs

* alignment file
* reference sequence

\vskip0.5em

Outputs

* VCF or BCF file

\vskip1em

Example\footnotesize

\hskip2em\texttt{bcftools mpileup -f ref.fa aln.bam | bcftools call -mv}

\vskip2em\normalsize
Tips\footnotesize

\vskip0.5em\hskip2em\texttt{bcftools mpileup}

\vskip0.1em\hskip3em - increase/decrease the required number (\texttt{-m}) and the fraction (\texttt{-F})
      of supporting reads for indel calling
\vskip0.1em\hskip3em - the \texttt{-Q} option controls the minimum required base quality (30)
\vskip0.1em\hskip3em - BAQ realignment is applied by default and can be disabled with \texttt{-B}
\vskip0.1em\hskip3em - streaming the uncompressed binary BCF (\texttt{-Ou}) is much faster than the default text VCF

\vskip0.5em\hskip2em\texttt{bcftools call}

\vskip0.1em\hskip3em - decrease/increase the prior probability (\texttt{-P}) to decrease/increase sensitivity

\vskip1em

\normalsize
General advice

* take time to understand the options
* play with the parameters, see how the calls change




# Factors to consider in calling

Many calls are not real, a filtering step is necessary

\vskip1em

False calls can have many causes

* contamination

* PCR errors

* sequencing errors
    * homopolymer runs

* mapping errors
    * repetitive sequence
    * structural variation

* alignment errors
    * false SNPs in proximity of indels
    * ambiguous indel alignment


# Callable genome

% ~/git/wtxt/logs/sandbox/other/ChangeLog   #1475659744

\small 

Large parts of the genome are still inaccessible

* the Genome in a Bottle high-confidence regions:
    * cover 89% of the reference genome
    * are short intervals scattered across the genome


\vskip1em

\centerline{\includegraphics[width=0.8\textwidth]{img/callable-chr1.pdf}}
\centerline{
    \includegraphics[width=0.4\textwidth]{img/callable-size.pdf}
    \includegraphics[width=0.4\textwidth]{img/callable-dist.pdf}
}

\begin{textblock}{8}(5,7.9)
\vonly[2-]{{\color{nRed}Q: Can you explain the spike at $\sim$100bp?}}
\vskip0.1em
\vonly[3]{{\color{nGreen}A: Short read sequencing}}
\end{textblock}


# Maximum depth

\centerline{\includegraphics[width=\textwidth]{img/artefact-depth.png}}

\vskip2em
\centering
{\color{nRed}Q: Why is the sequencing depth thousandfold the average in some regions?}
\vskip0.1em
\pause
{\color{nGreen}A: The reference genome is not complete. This
    sample was sequenced to 30x coverage, we can infer it has 
    $\sim$30 copies of this region.}




# Mapping errors

\centerline{\includegraphics[width=\textwidth]{img/aln-artefact-paralog.png}}

\vskip3em

\noindent{\color{nRed}Q: RNA-seq (top) and DNA data (bottom) from the same sample has been mapped onto the reference genome. Can you explain the novel SNVs?}
\pause
\vskip0.1em
\noindent{\color{nGreen}A: The reads originate in a paralog with 92\% identity.}



# Strand bias

\vonly[1]{\centerline{\includegraphics[width=\textwidth]{img/strand-bias-1.jpg}}}
\vonly[2]{\centerline{\includegraphics[width=\textwidth]{img/strand-bias-2.jpg}}}

\vskip3em

\vonly[1-]{{\color{nRed}Q: Is this a valid call?}\vskip0.1em}
\vonly[1]{\phantom{Q:}}
\vonly[2]{{\color{nGreen}A: No, it is a mapping artefact, the call is supported by forward reads only.}}


# Variant distance bias

\vonly[1]{\centerline{\includegraphics[width=\textwidth]{img/artefact-splice-site-vdb-1.png}}}
\vonly[2]{\centerline{\includegraphics[width=\textwidth]{img/artefact-splice-site-vdb-2.png}}}

\vskip2em
\centerline{\color{nRed}Q: Can you explain what happened here?}
\pause
\centerline{\color{nGreen}A: Processed transcript with introns spliced out.}


# Reproducibility

\centerline{\includegraphics[width=\textwidth]{img/artefact-reproducibility.png}}



# False SNPs caused by incorrect alignment

Pairwise alignemnt artefacts can lead to false SNPs 

* multiple sequence alignment is better, but very expensive

* instead: base alignment quality (BAQ) to lower quality of misaligned bases

\vskip3em
\vonly[1]{\centerline{\includegraphics[width=0.85\textwidth]{img/pileup-baq1.pdf}}}
\vonly[2]{\centerline{\includegraphics[width=0.85\textwidth]{img/pileup-baq2.pdf}}}

\vskip3em
\centerline{\color{nRed}Q: How many SNPs are real?}
\pause
\centerline{\color{nGreen}A: None.}





# How to estimate the quality of called SNPs?

Transitions vs transversions ratio, known as ts/tv

* transitions are 2-3$\times$ more likely than transversions

\vskip2em

\hskip2em\includegraphics[width=4.2\TPHorizModule]{img/Transitions-transversions-v4.pdf}



# Detour: Some causes of SNPs

Spontaneous chemical processes which lead to base modification or loss

* Deamination
    * methylated CpG dinucleotides: 5-methylcytosine $\rightarrow$ T
    <!-- methylation changes gene expression, methylcytosines can spontaneously deaminate to thymine, 4x less frequent than by chance, most common single nucleotide mutation -->

    \vskip0.5em
    * hydrolytic deamination of C $\rightarrow$ U \footnotesize (400 cytosines daily in each cell)\normalsize

    \vskip0.5em
    * A $\rightarrow$ hypoxanthine (pairs with C, A-to-G mutation)

\vskip1em

% See Lindahl above 
* Depurination (loss of A or G)
    * purines are cleaved from the backbone \footnotesize ($10^2$-$10^3$ daily in each cell)\normalsize
    * if base excision repair fails, random base is inserted

\vskip1em
DNA damage by mutagens

* base analogs 
    * incorporation of chemicals with different properties

* base-modifying agents

\vskip1em
Radiation

\begin{textblock}{1}(7.5,1.2)
\includegraphics[width=1.8\TPHorizModule]{img/cytosine-to-thymine.jpg}
\end{textblock}

\begin{textblock}{1}(7.5,2.4)
\includegraphics[width=1.5\TPHorizModule]{img/cytosine-to-uracil.jpg}
\end{textblock}

\begin{textblock}{1}(7.5,3.5)
\includegraphics[width=1.8\TPHorizModule]{img/adenine-to-hypoxanthine.pdf}
\end{textblock}

\begin{textblock}{1}(7.7,4.5)
\centerline{\includegraphics[width=2.8\TPHorizModule]{img/depurination.jpg}}
\end{textblock}

%  Wikipedia is probably wrong? 
%   G -> X -> A     https://en.wikipedia.org/wiki/Deamination
%   G -> X -> G     http://biowiki.ucdavis.edu/Genetics/Unit_II%3A_Replication,_Maintenance_and_Alteration_of_the_Genetic_Material/7._Mutation_and_repair_of_DNA/7.2_Reaction_with_mutagens
%   Both hypoxanthine and xanthine preferentially base-pair with cytosine, so G->X is not mutagenic
%   Lindahl, T 1993, Instability and decay of the primary structure of DNA, doi:10.1038/362709a0 



# Some causes of MNPs

UV-induced mutations (CC $\rightarrow$ TT in skin cells)

\vskip1em

\centerline{\includegraphics[width=\textwidth]{img/uv-mutation-igv-bezi3-chr12.jpg}}
\begin{textblock}{1}(7.5,4)
\includegraphics[width=2\TPHorizModule]{img/evil-sun.png}
\end{textblock}



# How to estimate the quality of called SNPs?

Transitions vs transversions ratio, known as ts/tv

* transitions are 2-3$\times$ more likely than transversions

\vskip2em

\hskip2em\includegraphics[width=4.2\TPHorizModule]{img/Transitions-transversions-v4.pdf}


\begin{textblock}{1}(6.3,0.09)
\includegraphics[width=2.9\TPHorizModule]{img/tstv-random.pdf}  % python tstv-normal-ancient.py
\vskip-0.6em
\includegraphics[width=2.9\TPHorizModule]{img/tstv-normal.pdf}  % python tstv-normal-ancient.py
\vskip-0.6em
\includegraphics[width=2.9\TPHorizModule]{img/tstv-ancient.pdf}  % python tstv-normal-ancient.py
\end{textblock}





# Indel calling challenges

The sequencing error rate is elevated in microsatellites

\vskip0.5em

Low reproducibility across callers

* 37.1% agreement between HapCaller, SOAPindel and Scalpel\newline \footnotesize Narzisi et al. (2014) Nat Methods, 11(10):1033

\vskip0.5em

Reads with indels are more difficult to map and align

* the aligner can prefer multiple mismatches rather than a gap
* indel representation can be ambiguous

\vskip1em 

% egypt.10AJ137.bam, 1:81688526-81688655 
\centerline{\includegraphics[width=10\TPHorizModule]{img/indel-realn-is-needed-1.png}}

\vskip1em

\vonly[1]{\centerline{\includegraphics[width=7\TPHorizModule]{img/indel-realn-is-needed-2.pdf}}}
\vonly[2]{\centerline{\includegraphics[width=7\TPHorizModule]{img/indel-realn-is-needed-3.pdf}}}



%   # Some causes of indels
%   
%   Replication slippage in short tandem repeats
%   
%   \vskip1em
%   \centerline{\hskip2em\includegraphics[width=6\TPHorizModule]{img/polymerase-slippage.pdf}}
%   
%   <!-- HD, polymerase slippage, hairpins: http://web.stanford.edu/group/hopes/cgi-bin/hopes_test/all-about-mutations/ -->
%   
%   \vskip1em
%   
%   Mutagents
%   
%   * intercalating agents




# Future of variant calling

Current approaches

* rely heavily on the supplied alignment, but aligners see one read at a time
* largely site based, do not examine local haplotype and linked sites

\vskip0.5em

Local *de novo* assembly based variant callers

* call SNPs, indels, MNPs and small SV simultaneously
* can remove alignment artefacts
* eg GATK haplotype caller, Scalpel, Octopus

\vskip0.5em

Variation graphs

* align to a graph rather than a linear sequence

\centerline{\includegraphics[width=0.5\textwidth]{img/colored-debruijn-graph.png}}
\centerline{\footnotesize Iqbal et al. (2012) Nat Gen 44(2):226}


# Single vs multi-sample and gVCF calling

VCF files can be **very** big, therefore we often store only variant 
sites\footnote{Annotated VCF with 3,781 samples, variant sites only, UK10k project $\ldots$ 680GB}


* however, variant-only VCFs are difficult to compare - was a site dropped because of a reference call
or because of low coverage?
* we need evidence for both variant and non-variant positions in the genome

\vskip0.5em

gVCF

* represents blocks of reference-only calls in a single record using the END tag
* symbolic allele in raw "callable" gVCFs allows to calculate genotype likelihoods
  only once (an expensive step), then do calling repeatedly as more samples come in

\centerline{\includegraphics[width=1\textwidth]{img/gvcf.pdf}}


# Functional annotation

VCF can store arbitrary INFO tags (per site) and FORMAT tags (per sample)

* describe genomic context of the variant (e.g. coding, intronic, UTR)
* predict functional consequence (e.g. synonymous, missense, start lost)

\vskip1em

Several tools for annotating a VCF, only few are haplotype-aware

\vskip0.1em\ibox[10em]{BCFtools/csq} \footnotesize http://github.com/samtools/bcftools \normalsize
\vskip0.1em\ibox[10em]{VEP Haplosaurus} \footnotesize http://github.com/willmclaren/ensembl-vep \normalsize

\vskip1em
\centerline{\includegraphics[width=0.8\textwidth]{img/compound-variants.pdf}}

