# Bioinformatics: Introduction and Methods

生物信息学：导论与方法

https://www.coursera.org/course/pkubioinfo

## Introduction and History

The miracle of life.

Genome: the "manual of life"

mitochondrial DNA, epigenetics, environments/nurture, chance

Human genome has 3.1 billion base pairs. ~2.9% of the bases encode genes. ~97% of the genome was previously called "junk". They contain the regulatory elements that encode instructions on when, where, which, and how much proteins to make.

[The Sequence Read Archive (SRA): Introduction](https://www.ncbi.nlm.nih.gov/sra/docs/)

The SRA is NIH's primary archive of high-throughput sequencing data and is part of the International Nucleotide Sequence Database Collaboration (INSDC) that includes at the NCBI Sequence Read Archive (SRA), the European Bioinformatics Institute (EBI), and the DNA Database of Japan (DDBJ). Data submitted to any of the three organizations are shared among them.

Opportunities and challenges hand‐in‐hand:
the driving forces of bioinformatics
- High‐throughput data
  - Huge amount
  - Explosive growth
  - Low signal‐to‐noise ratio
  - Multiple types
- Requirements for the methods
  - Data needs to be stored in efficient **ontology**‐based **database** systems
  - The huge amount of data requires **efficient** algorithms
  - Exponential growth requires **scalable** methods
  - The low signal‐to‐noise ratio requires **accurate** methods
  - Multiple types of data require data **integrative** methods

What is bioinformatics?

Life sciences -- Computer/Computational sciences

### Bioinformatics defination

Bioinformatics: an interdisciplinary field that develops and applies computer and computational technologies to study biomedical questions.
- As a technology, bioinformatics is a powerful technology to manage, query, and analyze big data in life sciences
- As a methodology, bioinformatics is a top-down, holistic, data-driven, genome-wide, and systems approach that generates new hypotheses, finds new patterns, and discovers new functional elements

![interdisciplinary](./figures/interdisciplinary.png)

### History of Bioinformatics

**The Bio- in Bioinformatics**

**Genotype ----> Phenotype**

DNA/Genome --> RNA --> Proteins --> Molecular Networks --> Cells --> Physiology/Disease

**Data ----> Discovery**

Data management --> Data computation --> Data mining --> Modeling/Simulation

History
- 1945, first electronic computer, ENIAC
- 1948, Shannon information theory
- 1953, **DNA double-helix structure**
- 1953, **first protein sequence determined**
- 1953, **first protein structure solved**
- 1953, first commerical (mass-produced) computers
- 1953, theory of games, grammars, FORTRAN
- 1962-1965, **molecular evolutionary clock**
- 1962-1965, cellular automata, theory of computation
- 1969, ARPANET created
- 1969, Email, Ethernet and TCP described
- 1977, **Sanger DNA sequencing**
- 1980-1990, personal computer, DNS launched
- 1990, **human genome project started**
- 1995, **H. influenzae genome sequenced**
- **E.coli, yeast, worm, fly genomes**
- 2001/2004, **human genomes sequenced**
- 2005, **Next-generation sequencing**
- 1951-1952, **1st computer program to determine protein structure**
- 1962, COMPROTEIN, M.O. Dayhoff
- 1965, protein atlas
- 1966, evolution of proteins
- 1967, phylogenetic tree
- 1970, first apperances of the word 'Bioinformatics': Hesper B., Hogeweg P. (1970). "Bioinformatica: een wekconcept", Kameleon, 1(6): 28-29. "The study of informatic processes in biotic systems"
- 1970, sequence alignment by DP
- 1971, **Protein Data Bank (PDB)**
- 1974, **Chou-Fasman 2nd structure prediction**
- 1975, **Levitt-Warshel computer simulation of protein folding**
- 1978, PAM substitution matrix, BLOSUM (1992)
- 1981, Smith-Waterman local alignment
- 1982, Genbank
- 1983, Protein Information Resource (PIR)
- 1990, 1997, BLAST, Gapped BLAST, PSI BLAST
- 1994, **CASP structure prediction assessment**
- Microarray gene expression analysis
- Gene prediction from genome
- Genome alignment
- 2002, BLAT
- 2007, SRA
- NGS data analysis

Examples of contributions to computer sciences
- Neural networks
- Genetic programming
- WAN

Current Bioinformatics Journals
- Bioinformatics
- BMC Bioinformatics
- BMC Systems Biology
- Briefings in Bioinformatics
- Bulletin of Mathematical Biology
- Cancer Informatics
- Computational Biology and Chemistry
- Computers in Biology and Medicine
- Database: The Journal of Biological Databases and Curation
- IEEE/ACM Transactions on Computational Biology and Bioinformatics
- Journal of Bioinformatics and Computational Biology
- Journal of Biomedical Informatics
- Journal of Computational Biology
- Journal of Integrative Bioinformatics
- Journal of Mathematical Biology
- Journal of Theoretical Biology
- PLoS Computational Biology
- Source Code for Biology and Medicine
- Statistical Applications in Genetics and Molecular Biology
Others:
- Nucleic Acids Research
- Genome Research
- Nature Methods
- Nature Biotechnology

### Bioinformatics in Mainland China

Rapid growth since 1990s. Driving force: 1. Internet 2. Genomics 3. increased research funding (Ministry of Science and Technology (MOST), Natural Science Foundation of China (NSFC)) 4. critical mass of researchers & students.

### Reading materials:

Required:

Bioinformatics. 2003 Nov 22;19(17):2176-90. Early bioinformatics: the birth of a discipline--a personal view. Ouzounis CA, Valencia A. PMID: 14630646
PLoS Comput Biol. 2011 Mar;7(3):e1002021. doi: 10.1371/journal.pcbi.1002021. Epub 2011 Mar 31. The roots of bioinformatics in theoretical biology. Hogeweg P. PMID: 21483479
PLoS Comput Biol. 2008 Apr 25;4(4):e1000020. doi: 10.1371/journal.pcbi.1000020. Bioinformatics in China: a personal perspective. Wei L, Yu J. PMID: 18437216
Optional:
Nature. 1991 Jan; 349(10):99 Towards a paradigm shift in biology. Gillbert W.

## Sequance Alignment

- Biology
  - What is the biological question or problem?
- Data
  - What is the input data?
  - What other supportive data can be used?
- Model
  - How is the problem formulated computationally?
  - Or, what’s the data model?
- Algorithm
  - What is the computational algorithm?
  - How about its performance/limitation?

### Essential Concepts

Biological Question: “How can we determine the similarity
between two sequences?”

Why is it important?
- Similar sequence -> Similar structure -> Similar function (The “Sequence‐to‐Structure‐to‐Function Paradigm”)
- Similar sequence -> Common ancestor (“Homology”)

Sequence Alignment in Biology

The purpose of a sequence alignment is to line up all residues in the inputted sequence(s) for **maximal level of similarity**, in the sense of **their functional or evolutionary relationship**.

[Pairwise Sequence Alignment](https://www.ebi.ac.uk/Tools/psa/)

[Sequence_alignment](Sequence_alignment.png)
- `|`: same residue
- `:`: denotes the similarity of different residue. It's similar.
- `.`: denotes the similarity of different residue. Not similar.
- `-`: gaps, insertion or deletion (indel), gap penalty
  - Affine gap penalty: **opening** a gap receives a penalty of d; **extending** a gap receives a penalty of e. So teh total Penalty for a gap with length n would be:
  Penalty = d + (n-1)*e
  - Final Score = (sum of substitution scores) + (-1)*(sum of Gap Penalty)

Substitution matrix: EBLOSUM62
- Symmetry
- Context-insensitive

### Readings

Required

J Mol Biol. 1981 Mar 25;147(1):195-7. Identification of common molecular subsequences. Smith TF, Waterman MS. PMID: 7265238(Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper. However you might also find other websites holding the full-text for this paper, because it is so famous.)

J Mol Biol. 1970 Mar;48(3):443-53. A general method applicable to the search for similarities in the amino acid sequence of two proteins. Needleman SB, Wunsch CD. PMID: 5420325 (Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper. However you might also find other websites holding the full-text for this paper, because it is so famous.)

Optional

Proc Natl Acad Sci U S A. 1992 Nov 15; 89(22): 10915–9. Amino acid substitution matrices from protein blocks. S Henikoff and J G Henikoff. PMID: 1438297

Nature. 1953 Apr 25;171(4356):737-8.Molecular structure of nucleic acids; a structure for deoxyribose nucleic acid. WATSON JD, CRICK FH. PMID: 13054692 (Note: this paper is so old, which seems to be the reason why PubMed does not provide a full-text link for this paper. Nevertheless, this paper is so famous, and you can find it easily on the Internet (and even on Wikipedia! Try it.) . Here is the link to the full-text of this paper on the website of Nature Publishing Group: http://www.nature.com/nature/dna50/watsoncrick.pdf)

Horizons in Biochemistry. 1962 pp. 189-225. Molecular disease, evolution, and genetic heterogeneity. Zuckerkandl E., Pauling L.B. (Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper. The format of reference here has not been confirmed to be correct yet, and we will check it later.)

Acta Chem. Scand. 1963;17:S9-S16. Pauling L., Zuckerkandl E. Chemical paleogenetics: molecular "restoration studies" of extinct forms of life (Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper.)

Evolving Genes and Proteins. 1965 pp. 97-166. Zuckerkandl E., Pauling L. Evolutionary divergence and convergence in proteins (Note: this paper is not archived in PubMed. You can, however, access the full-text of this paper from here: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.550.4751&rep=rep1&type=pdf

Kameleon. 1970 1(6);28-29. Hesper B., Hogeweg P. Bioinformatica: een werkconcept.

Simulation. 1978 Sep 1;31(3):90. Hogeweg, P. Simulating the growth of cellular forms. (Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper.)

Acta Cryst. 1952;5:109-116. Bennett J.M., Kendrew J.C. The computation of fourier syntheses with a digital electronic calculating machine.(Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper.)

Proceedings-Fall Joint Computer Conference. 1962. Dayhoff M., Ledley R. COMPROTEIN: a computer program to aid primary protein structure determination. (Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper.)

Science. 1966 Apr 15;152(3720):363-6. Eck RV, Dayhoff MO. Evolution of the structure of ferredoxin based on living relics of primitive amino Acid sequences. PMID: 17775169 (Note: you might need a subscription to the Science magazine in order to view the full-text of this paper.)

Science. 1967 Jan 20;155(3760):279-84. Fitch WM, Margoliash E.Construction of phylogenetic trees. PMID: 5334057 (Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper.)

J Mol Biol. 1970 Mar;48(3):443-53. A general method applicable to the search for similarities in the amino acid sequence of two proteins. Needleman SB, Wunsch CD. PMID: 5420325 (Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper.)

Atlas of protein sequence and structure. 1978. Chapter 22: A model of evolutionary change in proteins. Dayhoff M. O., Schwartz R. M. (Note: this paper is not archived in PubMed. I (Yang Ding the TA) cannot find the full-text of paper online as I'm writing this page, and since there are copyright issues, I recommend you to refer to some library to find the full-text of this paper.)

J Mol Biol. 1981 Mar 25;147(1):195-7. Identification of common molecular subsequences. Smith TF, Waterman MS. PMID: 7265238 (Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper. However you might also find other websites holding the full-text for this paper, because it is so famous.)

J Mol Biol. 1990 Oct 5;215(3):403-10. Basic local alignment search tool. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. PMID: 2231712 (Note: you might need a subscription to the Journal of Molecular Biology magazine in order to view the full-text of this paper. However you might also find other websites holding the full-text for this paper, because it is so famous.)

Nucleic Acids Res. 1997 Sep 1;25(17):3389-402. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman DJ. PMID: 9254694

Bioinformatics. 2003 Nov 22;19(17):2176-90. Early bioinformatics: the birth of a discipline--a personal view. Ouzounis CA, Valencia A. PMID: 14630646

## Global Alignment by Dynamic Programming

Pairwise Sequence Alignment: in Maths
- Input data:
  - Two sequences S1 and S2
- Parameter(s)
  - A scoring function f for
    - Substitutions
    - Gaps
- Output:
  - The optimal alignment of S1 and S2, which has the maximal score.
  $\arg \max_{ali}(f(ali(S1, S2))$

Sequence Alignment: Enumerate?

$\frac{(2n)!}{(n!)^2}$

Sequence Alignment:

What is the computational Algorithm? A residue can either
- Align to other residue, or
- Align to a gap

The **best alignment** that ends at a given pair of symbols is the **best alignment** of the sequences up to that point, plus the **best alignment** for the two additional symbols.

Dynamic Programming

Dynamic Programming solves problems by combining the solutions to sub‐problems.
- Break the problem into smaller sub‐problems
- Solve these sub‐problems optimally recursively
- Use these optimal solutions to construct an optimal solution for the original problem

The **best alignment** that ends at a given pair of symbols is the **best alignment** of the sequences up to that point, plus the **best alignment** for the two additional symbols.

![Sequence_alignment_dp](./figures/Sequence_alignment_dp.png)

![Sequence_alignment_dp1](./figures/Sequence_alignment_dp1.png)

![Sequence_alignment_dp2](./figures/Sequence_alignment_dp2.png)

![Sequence_alignment_dp3](./figures/Sequence_alignment_dp3.png)

![Sequence_alignment_dp4](./figures/Sequence_alignment_dp4.png)

![Sequence_alignment_dp5](./figures/Sequence_alignment_dp5.png)

![Sequence_alignment_dp6](./figures/Sequence_alignment_dp6.png)

![Sequence_alignment_dp7](./figures/Sequence_alignment_dp7.png)

![Sequence_alignment_dp8](./figures/Sequence_alignment_dp8.png)

![Sequence_alignment_dp9](./figures/Sequence_alignment_dp9.png)

![Sequence_alignment_dp10](./figures/Sequence_alignment_dp10.png)

![Sequence_alignment_dp11](./figures/Sequence_alignment_dp11.png)

![Sequence_alignment_dp12](./figures/Sequence_alignment_dp12.png)

Homology
- derived from a common ancestor
  - ortholog: derived from speciation
  - paralog: derived from duplication
- Ortholog comes with speciation. Paralog comes with duplication

![similarity_vs_identity](./figures/similarity_vs_identity.png)

Similarity Matrix
- For nucleotides,
  - usually only distinguish
match / mismatch
(identity matrix)
for sequence alignment
  - but a more complicated
substitution model is used
for phylogeny reconstruction
- For amino acids,
  - PAM (1978, Margaret Dayhoff)
    - Two sequences are 1 PAM apart
    if they differ in 1 % of the residues.
    - 1 PAM = one step of evolution
- BLOSUM (1992, Steven Henikoff & Jorja Henikoff)
  - computed by looking at "blocks" of conserved sequences found in multiple protein alignments

How to let computer do this job?
- How to measure similarity?
  - Similarity matrix
- How to find out alignment?
  - Dot matrix
  - Dynamic programming
  - BLAST

## From Global to Local

Global Alignment: End‐to‐end

Needleman–Wunsch algorithm

Identify similar sub‐sequence

![global_to_local](./figures/global_to_local.png)

![global_to_local1](./figures/global_to_local1.png)

![global_to_local2](./figures/global_to_local2.png)

![global_to_local3](./figures/global_to_local3.png)

![global_to_local4](./figures/global_to_local4.png)

![global_to_local5](./figures/global_to_local5.png)

![global_to_local6](./figures/global_to_local6.png)

![global_to_local7](./figures/global_to_local7.png)

![global_to_local8](./figures/global_to_local8.png)

## Alignment with Affine Gap Penalty

![global_to_local9](./figures/global_to_local9.png)

![global_to_local10](./figures/global_to_local10.png)

![global_to_local11](./figures/global_to_local11.png)

![global_to_local12](./figures/global_to_local12.png)

![global_to_local13](./figures/global_to_local13.png)

![global_to_local14](./figures/global_to_local14.png)

Smith-Waterman Algorithm
- Compared with global alignment：
  - Zero could terminate the current local alignment
  - Mismatch must be negative scored
- Other properties：
  - Suitable to identify conserved local sequence (substring)
  - Guaranteed to find the best local alignment
  - Perform poorly when dealing with separated regions within a long sequence

## Sequence Database Search

Sequence Database Searching
- Rather than do the alignment pair‐wise, it’s more often to search **sequence database** in a **high‐throughput** style.
- Or, identify similarities between
  - **novel query sequence**
  whose structures and functions are usually unknown and/or uncharacterized
  - **sequences in (public) databases**
  whose structures and functions have been elucidated and annotated.
- The **query sequence** is compared/aligned with every sequence in the database
- **Statistically significant hits** are assumed to be related to the query sequence
  - Similar **function/structure**
  - Common **evolutionary ancestor**

BLAST: Intro
- To make the alignment effectively, a Heuristic algorithm BLAST
(**B**asic **L**ocal **A**lignment **S**earch **T**ool) is proposed by Altschul et al in 1990.
- BLAST finds the highest scoring locally optimal alignments between a query sequence and a database.
  - Very **fast** algorithm
  - Can be used to search **extremely large** databases
  - Sufficiently **sensitive** and **selective** for most purposes
  - **Robust** – the default parameters just work for most cases

[UniProt](https://www.uniprot.org/)

## BLAST Algorithm: a Primer

BLAST Ideas: Seeding‐and‐extending
1. Find matches (**seed**) between the query and subject
2. Extend seed into High Scoring Segment Pairs (**HSP**s)
  Run Smith‐Waterman algorithm on the specified region only.
3. Assess the reliability of the alignment.

Seeding: For a given word length w (usually 3 for proteins and 11 for nucleotides), slicing the query sequence into multiple continuous “seed words”.

Speedup: Index database
The database was pre‐indexed to quickly locate all positions in the database for a given seed.

Speedup: mask low‐complexity
Low complexity sequences yield false positives.

To improve sensitivity, in addition to the seed word itself, the BLAST also use these highly similar “neighbourhood words” (based on the substitution matrix) for seeding.

Quality Assessment
Given the large data volume, it’s critical to provide some measures for assessing the **statistical significance** of a given hit.

E‐Value: How a match is likely to arise by chance
- The expected number of alignments with a given score that would be expected to occur at random in the database that has been searched
  - e.g. if E=10, 10 matches with scores this high are expected to be found by chance
  $E = kmne^{-\lambda S}$
  $p = 1-e^{-E}$

  m: length of the query sequence
  n: length of the database
  S: score of HSP
  E>1, the alignment occurred by chance
  E<0.1 or 0.05, statistically significant
  E<10-5, high similarity

Heuristic (pronounced hyu‐RIS‐tik, Greek: "Εὑρίσκω", "find" or "discover") refers to experience‐based techniques for problem solving, learning, and discovery. (Source: Wikipedia)

not best, but good enough.

- Key heuristics in BLAST
  - Seeding‐and‐extending: looking for seeds of high scoring alignments ONLY
  Use dynamic programming selectively
- Tradeoff: speed vs. sensitivity
  - Empirically, 1000 ~ 10000 times faster than plain Dynamic‐Programming‐based local alignment
  - But suffer from low sensitivity, especially for distant
sequences (e.g. E.coli -> human)

Why BLAST?

"**Homology** is the central concept for all of biology." -- David Wake. Science, 1994

BLAST is the tool most frequently used for calculating sequence **similarity**, by searching the databases.

If you work with one or a few proteins or genes, it can tell you about their conservation, active sites, structure and regulation in other organisms, etc.

What BLAST does?
- Identity: the occurrence of exactly the same nucleotide or amino acid in the same position in aligned sequences.
- Similarity: measure the sameness or difference of the sequences 
- Homology: is defined in terms of shared ancestors. Homologous sequences are often similar. Sequence regions that are homologous are also called **conserved** regions.

![alignment_algorithm](./figures/alignment_algorithm.png)

![blast](./figures/blast.png)

![blast1](./figures/blast1.png)

### Readings

Required

J Mol Biol. 1990 Oct 5;215(3):403-10. Basic local alignment search tool. Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. 

(Note: currently there is no full-text link for this article in PubMed. However you might also find other websites holding the full-text for this paper, because it is so famous.)

http://www.fi.muni.cz/~lexa/blast.html

Nucleic Acids Res. 1997 Sep 1;25(17):3389-402. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Altschul SF, Madden TL, Schäffer AA, Zhang J, Zhang Z, Miller W, Lipman, D.J. PMID: 9254694

Genome Biol. 2001;2(10):REVIEWS2002. Having a BLAST with bioinformatics (and avoiding BLASTphemy). Pertsemlidis A, Fondon JW 3rd. PMID: 11597340

Nucleic Acids Res. 2006 Jul 1;34(Web Server issue):W6-9. BLAST: improvements for better sequence analysis. Ye J, McGinnis S, Madden TL. PMID: 16845079

Optional

Nat Genet. 1993 Mar;3(3):266-72. Identification of protein coding regions by database similarity search. Gish W, States DJ. PMID: 8485583

Note: this paper is supplementary and thus OPTIONAL. You might need a subscription to the Nature Genetics magazine in order to view the full-text of this paper.

Proc Natl Acad Sci U S A. 1990 Mar;87(6):2264-8. Methods for assessing the statistical significance of molecular sequence features by using general scoring schemes. Karlin S, Altschul SF. PMID: 2315319

Note: this paper is supplementary and thus OPTIONAL.

Proc Natl Acad Sci U S A. 1993 Jun 15;90(12):5873-7. Applications and statistics for multiple high-scoring segments in molecular sequences. Karlin S, Altschul SF. PMID: 8390686

Note: this paper is supplementary and thus OPTIONAL.

Methods Enzymol. 1996;266:554-71. Analysis of compositionally biased regions in sequence databases. Wootton JC, Federhen S. PMID: 8743706

Note: this paper is supplementary and thus OPTIONAL. You might need a subscription to the Methods in Enzymology magazine in order to view the full-text of this paper.

Comput. Chem. 1993 Jun;17(2):149-163. Statistics of local complexity in amino acid sequences and sequence databases. Wootton JC, Federhen S. Full-text at ScienceDirect

Note: this paper is supplementary and thus OPTIONAL. You might need a subscription to the Computers & Chemistry magazine in order to view the full-text of this paper.

Recommended http://www.ncbi.nlm.nih.gov/pmc/articles/PMC290255/ (Recommended by Sebastien Beaune)

Readings

Required

Nature Biotechnology 22, 1315 - 1316 (2004). What is a hidden Markov model? Sean R Eddy. http://www.fing.edu.uy/~alopeza/biohpc/papers/hmm/Eddy-What-is-a-HMM-ATG4-preprint.pdf

Optional

Nucleic Acids Res. 1998 Jan 1;26(1):320-2. Pfam: multiple sequence alignments and HMM-profiles of protein domains. Sonnhammer EL, Eddy SR, Birney E, Bateman A, Durbin R. http://nar.oxfordjournals.org/content/26/1/320.full.pdf+html

## Markov Model

Andrei Andreyevich Markov (1856-1922). A Markov chaini describes a discrete stochastic process at successive times. The transitions from one state to any of all states, including itself, are governed by a probability distribution.
- $P(X_t | X_1...X_{t-1}) = P(X_i | X_{t-m}...X_{t-1})$
  - $X_t = f(X_{t-1}, X_{t-2}, ..., X_{t-m})$
- A chain of random variables in whick the next one depends (only) on the current one
  - $-P(X_t | X_{1}...X_{t-1}) = P(X_{i} | X_{i-1})$

Transition Probability
- $a_{kl} = P(X_{t} = S_{l} | X_{t-1} = S_{k})$, state k transite to state l.
- $a_{lk} = P(X_{t} = S_{k} | X_{t-1} = S_{l})$, state l transite to state k.
- $a_{kl}$ is not necassarily equal to $a_{lk}$

Hidden Markov Model (HMM)
- In addition to State Transition Probability, each state of HMM has a probability distribution over the possible output tokens (Emission Probability)
- Thus, a HMM is consist of two strings of information
  - The state path, which is not directly visible
  - The token path (emitted sqquence), based on the observable token path, we can infer the underling state path
- Transition and Emmision probability.
- Given a HMM, a sequence of tokens could be generated as following:
  - when we "visit" a state, we emit a token from the state's emission probability distribution
  - then, we choose which state to visit next, according to the state's transition probability distribution
    - transition probability, $a_{kl} = P(x_{t} = S_{l} | x_{t-1} = S_{k}$
    - emission probability, $e_{k}(b) = P(y_{i} = b | x_{i} = S_{k})$
  - $P(X, Y) = \prod_{i=1}^L(e_{x_{i}}(y_{i}*a_{x_{i}x_{i+1}})$
- M: Match (not necessarily identical); X: Insert at sequence X (delete at sequence Y); Y: Insert at sequence Y (delete at sequence X), $\delta$: Gap open; $\epsilon$: Gap extension.
  ||M|X|Y|
  |:------:|:------:|:------:|:------:|
  |**M**|$1-2\delta$|$\delta$|$\delta$|
  |**X**|$1-\epsilon$|$\epsilon$|0|
  |**Y**|$1-\epsilon$|0|0|

Sequence alignment with HMM
- Each "token" of the HMM is an aligned pair of two residues (M state), or of a residue and a gap (X or Y state)
  - transition and emission probabilityes define the probability of each aligned pair of sequences
- Based on the HMM, each alignment of two sequences can be assigned with a probability
  - given two input sequences, we look for an alignment with the maximum probability
    - arg max(P(S1, S2, ali))
    - $P(X, Y, ali) = max(P_{M}(n,m), P_{X}(n,m), P_{Y}(n,m))$

Given the nature of HMM, many different state paths can give rise to the same token sequence, so we can simply sum up them together to get the full probability of a given token sequence: $P(X, Y) = \sum_{ali}(X, Y, ali)$

HMM: as a predictor

The Most Simple Gene Predictor (MSGP): Given a stretch of genomic sequence, where are the coding regions and where are noncoding regions?

Training the model
- What we need to train?
  - transition probabilities between states
  - emission probabilities for each state
- Estimate probabilities from known "training set"
  - an annotated genomic region, with coding/noncoding sequences labeled
- $P_{coding}(i+1) = e_{coding}(x_{i+1}) max_{k \in (coding, noncoding)}(P_{k}(i)a_{k->coding})$
- $P_{noncoding}(i+1) = e_{noncoding}(x_{i+1}) max_{k \in (coding, noncoding)}(P_{k}(i)a_{k->noncoding})$
- $P(X, S) = max(P_{coding}(n), P_{noncoding}(n))$
- Logarithmic transformation: ease calculation: $Log(a*b) = Log(a) + Log(b)$
- "CGAAAAAATCG" will get the result "NNCCCCCCNNN" usign MSGP.

GenScan:
- Chris Burge (1996): A 27-state sem-HMM
- A simpler model: 19-state
- A model taking UTR introns into account: 35-state

By decoupling states and tokens, Hidden Markov Model (HMM) provides a sound probability framework to model complex biological sequences.

Readings Required

Nature Biotechnology 22, 1315 - 1316 (2004). What is a hidden Markov model? Sean R Eddy. http://www.fing.edu.uy/~alopeza/biohpc/papers/hmm/Eddy-What-is-a-HMM-ATG4-preprint.pdf

Optional

Nucleic Acids Res. 1998 Jan 1;26(1):320-2. Pfam: multiple sequence alignments and HMM-profiles of protein domains. Sonnhammer EL, Eddy SR, Birney E, Bateman A, Durbin R. http://nar.oxfordjournals.org/content/26/1/320.full.pdf+html

week 6

## Next Generation Sequencing (NGS)

Read: A short DNA fragment which is 
a readout by sequencer (in FASTQ format):
- DNA sequence (symbols)
- Quality information

Quality: Given p = the probability of a base calling is wrong, its Quality Score can be written as:
- Q = $-10*\log_{10}(p)$

Association Study: for the given phenotypic trait, "functional variants" could be identified by comparing allele frequencies at hundreds of thousands of polymorphic sites, i.e allele A is associated with phenotypic trait P if (and only if) people who have P also have A more (or less) often than would be predicted from individual frequencies of A and P in the assessed population.

RNA-Seq: Explore the transcriptome. "A transcriptome is a collection of all the transcripts present in a given cell." (NHGRI factsheet, NIH, US)

Chromatin ImmunoPrecipitation Sequencing (ChIP-Seq): Profile Protein-DNA interaction. (deep sequencing technology)

Reads Mapping and Variants Calling: the first step of deep sequencing

Mapping: Input data
- Reference Genome
  - Nucleotide
  - Length: humdreds of Mb per chromosome
  - ~3 Gb in total (for human genome)
- Reads
  - Nucleotide, with various qualities (relatively high error rate: 1e-2~1e-5)
  - Length: 36~80 bp per read
  - Hundreds of Gbs per run

 "Embeded" Alignment: One sequence is "embeded" in the other sequence (NGS Reads, PCR primer, etc.) What we need here is actually a hybrid "global-local" alignment.
 - "Global" for short sequence (i.e. NGS Read)
 - But "Local" fro long sequence (i.e. Reference Genome)
 - In particular, the surrounding "overhang" gaps should be not penalized

BLAST Ideas: Seeding-and-extending
1. Find matches (seed) between the query and subject 
2. Extend seed into High Scoring Segment Pairs (HSPs) 
   - Run Smith-Waterman algorithm on the specified region only
3. Assess the reliability of the alignment

The speed is a big mater. Using Prefix tree or Suffix tree. 

HBS: A native hash function

Pigeonhole principle

Mapping Quality: Given reference sequence z (length L), a read sequence x (length l), u is the alignment position of x on z, the probability that z actually coming from the position u is p(z|x,u)

$p(z|x,u) = \prod_{mismatch}p(z_{i})$

$SQ(U) = \log (p(z | x, u)) = \sum_{mismatch}p(z_{i}) = \sum_{mismatch}Q(z_{i})$

If we assume that a uniform NULL model, i.e. the read can randomly come from all possible positions with equal probability, then the error of mapping to a specified position u could be written as:

$E(u) = \frac{SQ(u)}{\sum_{i}SQ(i)}$

Genetic Variants
- SNV: single nucleotide variant
  - Substitution (SNP)
  - Indel: insertion/deletion
- Structural Variation (SV)
  - Large-scale insertion/deletion
  - Inversion
  - Translocation
  - Copy Number Variation (CNV)

SNP calling is not genotyping: source: Nature Reviews Genetics 12, 443-451.
- SNP calling aims to determine in which positions there are polymorphisms or in which positions at least one of the bases differs from a reference sequence
- Genotype calling is the process of determining the genotype for each individual and is typically only done for positions in which a SNP or a 'variant' has already been called

Counting: an intuitive (and naive) approach
- counting high-confident, non-reference allele (i.e. Quality >= 20)
  - Freq < 20% or > 80%: homozygous genotype
  - Otherwise: heterozygous
- Works well for "deeply sequenced regions" (DSR), i.e. depth > 25x
  - But suffer from under-calling of heterozygous genotypes for low-coverage regions
  - And can't give an objective measurement for reliability

  Source: Nature Reviews Genetics 12, 443-451.

  A simple probabilistic model for genotyping
  1. For a diploid genome, there will be at most two different alleles (A and a) observed at a given site:
    - 3 possible genotypes: <A,A>, <A,a>, <a,a>
    - number of A: k, number of a: n-k
2. Then, the probability for each genotypes is
  - P(D | <A,A>) = the probability that we have (n-k) sequencing errors at this site $\prod_{n-k}P(x_{i})$
  - Similarly, we can see the P(D | <a,a>) = $\prod_{k}P(x_{i})$
  - P(D | <A,a>) = 1- (P(D | <A,A>) + P(D | <a,a>))
3. Bayes Formula can be further employed to calculate posterior probabilities, i.e. P(<A,A> | D), P(<a,a> | D), P(<A,a> | D) if we can estimate the prior probabilities P(<A,A>), P(<a,a>) and P(<A,a>).

outline:
- BWA & BWT (Burrows-Wheeler transform) algorithm
  - BWT is the compression algorithm used in BWA 
  - Lossless compression
  - Sort and transform the char matrix with string rotation
  - Reverse-char method was utilized for match
  - Cannot handle gap
- Variant caller
  - samtools: mpileup + bcftools
  - GATK (Genome Analysis ToolKit): UnifiedGenotyper, HaplotypeCaller
  - [ref](http://www/broadinstitute.org/gatk/guide/best-practices)
- Demonstration

Introduction of Likelihood and Bayesian approach

Genotyper of MAQ and SNVMix

Readings Required

Nat Biotechnol. 2009 May;27(5):455-7. How to map billions of short reads onto genomes. Trapnell C, Salzberg SL. PMID: 19430453

Nat Rev Genet. 2010 Jan;11(1):31-46. Sequencing technologies - the next generation. Metzker ML. PMID: 19997069

An Introduction to Hidden Markov Models for Biological Sequences by Anders Krogh. Optional Nat Rev Genet. 2011 Jun;12(6):443-51. Genotype and SNP calling from next-generation sequencing data. Nielsen R, Paul JS, Albrechtsen A, Song YS. PMID: 21587300

Nat Rev Genet. 2013 Jul;14(7):460-70. Sequencing studies in human genetics: design and interpretation. Goldstein DB, Allen A, Keebler J, Margulies EH, Petrou S, Petrovski S, Sunyaev S. PMID: 23752795

## Variant Databases

Where did our genetic variations come from?
- inherited from parents
- de novo mutations
- somatic mutations

Types of genetic variations in a human genome
- Chromosomal aneuploidy
- Structural Variations (SVs): deletion, tandem duplication, interspersed duplication (Copy Number Variations (CNVs)), mobile-element insertion, novel sequence insertion, inversion, translocation.
- Indel - short Insertion/Deletion
  - Within intergenic/intronic regions
    - Non-frameshifting
    - Frameshifting
- SNV - Single nucleotide variation: there are about 3 million SNVs in one person's genome, equicalent to ~1/1000 frequency. promoter variant, coding variant, UTR variant, intronic variant, intergenic variants, non-coding RNA variant.
  - SNVs within coding regions:
    - Stop codon gain (nonsense)
    - Non-synonymous (nissense)
    - Synonymous (same sense/silent)
    - Affect splicing
    - Stop codon loss
  - Nonsense SNVs are usually considered deleterious
    - even though it is not always the case...
  - synonymous, intronic, and intergenic variations are often ignored
    - however, according to GWAS studies, 88% of trait-associated variants of weak effect are non-coding
    - they remain under-studied and new methods are needed
  - Most research so far had focused on missense SNVs
  - Known deleterious mutations are enriched in missense mutations
    - ~50 of all known mutations of Mendelian disorders are missense mutations. asscertainment bias?
  - Although missense SNVs change the protein sequence, many do not cause phenotypic changes.

Nomenclature
- Mutation: when the monor allele has a frequency less than 1% in the general population, we usually called it a mutation
- Polymorphism,: otherwise it is usually called a polymorphism. Sometimes the cut off of 5% is used.
- Variation/variant: mutations and polymorphisms together are called variations or variants.

Functional/Phenotypic "effects" of human genetic variations
- disease vs. normal
- deleterious vs. neutral
- personal trait differences (e.g., height)

- animal model phenotypic changes
- cellular phenotypic changes
 
 - protein function changes
 - protein structure changes

 Statistical and stochastic, not deterministic. Observations, not "truth"

Variant databases
- 1998 [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/): build 138 contains genetic variations from 131 species
  - SNPs
  - Indels
  - multinucleotide polymorphism
  - microsatellite markers
  - short tandem repeats
  - heterozygous sequences
- 2010 dbVar
- 2012 [1000 genomes](http://www.1000genomes.org/)
- 1987 [OMIM](http://www.omim.org/): Online Mendelian Inheritance in Man
- 1996 [HGMD](http://www.hgmd.cf.ac.uk/): Human Gene Mutation Database, a comprehensive collection of gene mutations that underlie, or are associated with, human genetic diseases, manually curated from literature.
- 2004 COSMIC
- 2007 LSDBs
- 2012 ClinVar

Conservation-based and rule-based methods: SIFT & PolyPhen
- 1999, an early attempt based on BLOSUM substitution matrix
- Mor successful methods since 2001
  - Conservation-based methods (e.g., SIFT) [Sort Intolerant From Tolerant substitutions](http://sift.jcvi.org/)
    - Important positions (such as active sites) tend to be conserved in the protein family across species.
      - mutations at well-conserved positions tend to be deleterious
    - Some positions have a high degree of diversity across species
      - mutations at these positions tend to be neutral
  - Rule-based methods (e.g., PolyPhen)
  - Machine learning classifier-based methods (e.g., [PolyPhen2](http://genetics.bwh.harvard.edu/pph2/), [SAPRED](http://www.sapred.cbi.pku.edu.cn/): **S**ingle **A**mino aicd **P**olymorphisms disease-association **Pred**ictor)

Sequence search database: SWISS-PORT; BLAST, FASTA, PSI-blast, CLUSTAL; MOTIF; SIFT & PolyPhen.

Homology Modeling: [SWISS-MODEL](http://swissmodel.expasy.org/)

Wherever there are challenges, there are opportunities.

Readings required

Bioinformatics. 2007 Jun 15;23(12):1444-50. Finding new structural and sequence attributes to predict possible disease association of single amino acid polymorphism (SAP). Ye ZQ, Zhao SQ, Gao G, Liu XQ, Langlois RE, Lu H, Wei L. PMID: 17384424

Optional:

Annu Rev Biophys Biomol Struct. 2000;29:291-325. Comparative protein structure modeling of genes and genomes. Martí-Renom MA, Stuart AC, Fiser A, Sánchez R, Melo F, Sali A. PMID: 10940251

week 7

## Mid-term tests

week 8

## Next generation sequencing: transcriptome analysis, and RNA-Seq

"A transcriptome is a collection of all the transcripts present in a given cell." (NHGRI factsheet, NIH, US)

Qualitative: identify all transcripts, i.e. expressed genes as well as their isoforms.

Quantitative: estimate the expression level of these transcripts, i.e. the transcript abundance of expressed genes/isoforms.

Real-Time qRT-PCR
- based on complementary hybridization
  - development of PCR technology
- widely accepted as the "Gold Standard"
- low-throughput
- Prior knowledge of transcript sequences needed!

Microarray
- DNA microarrays are used to analyze gene expression based on complementary hybridization reaction.
- labeled targets: RNAs derived from biological samples
- probes: a large number of ordered sets of immobilized nucleotide molecules with known sequences.
- prior knowledge of transcript sequences needed!

Expressed Sequence Tag (EST)
- Randomly sequencing individual clone of a cDNA library
  - short "tag": 60~500bp
  - one-shot: rather error-prone (~1%)
- NO prior knowledge of transcript sequences needed!
  - Not only measure, but also discovery

Random sampling the transcriptome: the detection power as well as sensitivity of RNA-Seq is highly dependent on the sequencing depth.
- 100~150x as a decent start for a typical mammalian transcriptome RNA-Seq.
- number of reads proportional to transcript abundance
- number of mapped reads proportional to transcript length
- number of mapped reads proportional to library depth

Genome DNA (Gene A vs Gene B (2:1))

Abundance of transcripts (2:2)

Observed reads count (4:2)

Observed abundance (4:2)

Genome DNA (Gene B vs Gene B (1:1))

Abundance of transcripts (2:2)

Total reads (depth): (2n:n)

Observed reads count (4:2)

Observed abundance (4:2)

From raw count to expression level

RPKM: the number of mapped Reads Per KB Per Million reads.

RPKM = $10^9\frac{C}{NL}$
- C: the number of mapped reads for specified transcript
- N: the number of total mapped reads
- L: the length fo the specified transcript

Mapping: Input Data

- Reference Genome
  - Nucleotide
  - Length: Hundreds of Mb per chromosome
  - ~3Gb in total (for human genome)
- Reads
  - Nucleotide, with various qualities (relatively high error rate: 1e-2~1e-5)
  - Length: 36~80 bp per read
  - Hundreds of Gbs per run

Handle junction reads: "Join exon" strategy
1. Build "conceptual junctions library" (CJL) for each known transcripts
2. Map RNA-Seq reads to the genome as well as the conceptual junction library
- Fast
- Can identify novel splicing isoforms 
- Can NOT find novel exons and novel genes

Handle junction reads: "Split reads" strategy
1. Unsplicingly map to genome, (Non-junction reads)
2. (if 1 failed), Split reads into several k-mer seeds
3. Retry with these seeds
4. Stitch mapped seeds together as whole read alignment 
- Slower than "join exon"
- Can identigy novel splicing isoforms
- Can find novel exons and novel genes

RNA-Seq: Mapping & Assembling
- Assembly: reconstruct full-length transcript sequences from the (mapped) reads
- Quantification: estimate the expression abundance for each transcripts

Assembling as a graph traveler

(Rationally) weighted edges based on various suppoeting evidences
- Existing transcript data: junction reads
- Sequence context: splicing junction sites, polyA signals, ...
- Existing annotation

Cufflinks: Transcript assembly, differential expression and differential regulation for RNA-Seq

Commandline tool:
- Tophat: Spliced reads mapper for RNA-Seq
- Cufflinks: Transcript assembly, differential expression, and differential regulation for RNA-Seq
- Cuffmerge: merge together several Cufflinks assemblies
- Cuffdiff
  - calculates expression in two or more samples: number of reads of each transcript proportional to its abundance
  - tests the statistical significance each observed chage in expression between them
R package:
- CummeRbund: An R package to aid and simplify the task of analyzing Cufflinks RNA-Seq output.
  - csDensity(genes(cuff_data))
  - csScatter(genes(cuff_data), 'C1', 'C2')
  - csVolcano(genes(cuff_data), 'C1', 'C2')
  - expressionBarplot(mygene)
  - expression Barplit(isoforms(mygene))

What can RNA-Seq do?
- Level 1
  - Gene expression (RPKM)
- Level 2
  - Transcriptome reconstruction & Alternative Splicing
  - Isoform abundance & Differential Expression
- Also
  - RNA-editing (DNA-SNP vs RNA-SNP)
  - eQTL(SNP-Expression)
  - Fusion

- Using result:
  - GO & Pathway Enrichment

Reading mapping, transcriptome reconstruction, differential expression.

Readings

Required

Nat Methods. 2011 Jun;8(6):469-77. Computational methods for transcriptome annotation and quantification using RNA-seq. Garber M, Grabherr MG, Guttman M, Trapnell C. PMID: 21623353

Nat Biotechnol. 2010 May;28(5):511-5. Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation. Trapnell C, Williams BA, Pertea G, Mortazavi A, Kwan G, van Baren MJ, Salzberg SL, Wold BJ, Pachter L. PMID: 20436464

week 9

## Transcriptome analysis with noncoding RNAs

### From information to knowledge

A non-coding RNA (ncRNA) is any RNA molecule that could function without being translated into a protein.

The DNA sequence from which a non-coding RNA is transcribed as the end product is often called an RNA gene or non-coding RNA gene.

Early discovered ncRNAs are mostly housekeeping
- "Assist"in translation in a necessary, but passive roles
- Constitutively expressed
- Include: rRNA, tRNA, snRNA, snoRNA, tmRNA, telomerase RNA...

Recently discovered regulatory ncRNAs since 2000
- actively regulate gene transcription and translation
- are involved in various gene regulations through multiple mechanisms
- many have specific expression patterns
- are widely encoded in the genome
  - the ENCODE (ENCyclopedia Of DNA Elements) pilot project suggested that over 90% of the human genome may be represented in primary transcripts.
  - over 95% of all the transcripts are noncoding. Some estimate the number of ncRNAs to be ~30000.

microRNA (miRNA)
- single-stranded RNAs of 21-23 (or some say 20-25) nt RNAs with regulatory functions when associated with a protein complex.
- in plants miRNAs can silence gene activity via destruction of homologous mRNA or blocking its translation. In animals, miRNAs inhibit translation by binding with imperfect homology to the 3' untranslated region of mRNA.

Xist: Beyond "small" ncRNA, X inactive-specific transcript

SCA8: Long ncRNA in human disease. SCA8 is mutated in one form of spinal cerebella ataxia.

Long ncRNAs
- Estimated ~2000+ in human
- Some, but not all, are mRNA-like, with Poly(A) tails.
- Most have unknown function. Many may function via cis or trans antisense pairing.
  - Dosage compensation (e.g. XIST)
  - Neuron development (e.g. SCA8)
  - Genetic imprinting (e.g. IGF/H19)
  - Post-transcriptional regulation
    - mRNA degradation or stabilization
  - Translational regulation
  - Modulate protein function by directly binding to the protein

### Data mining: identify long ncRNAs

Identification

Features ~ property of an entity
- structural features
- evolutionary features

Sequence features only
Mechanism neutral: works for both long and small ncRNAs
Accurate and fast

Initialized feature set
- properties of entity
- speculate based on existed knowledge
- certain statistic established by predecessors
- the data that is thought to be relevant

(Conceptually) Translated Product

Coverage: Coverage(S) = $\frac{L_{ORF} - (L_{mismatch} + 2 * L_{frameshift})}{Total\, Length}$

ORF integrity: indecates whether the predicted ORF begins with a start codon and ends with an in-frame stop condon

LOG-ODD score: indicator of the quality of a predicted ORF. The higher the score, the better the quality of the ORF.

Homologous

number of BLASTX hits: A true protein-coding transcript is likely to have more hits with known proteins than a non-coding transcript does

Hit Score: For a true protein-coding transcript, the hits are also likely to have higher quality 
$S_{i} = mean_{j}{-log_{10}{E_{ij}}}$, $i \in [0, 1, 2]$
HIT SCORE = $mean_{i \in {0, 1, 2}}{Si} = \frac{\sum_{i=0}^{2}S_{i}}{3}$

Frame Score: for a true protein-coding transcript, most of the hits are likely to reside within one frame, whereas for a true non-coding transcript, even if it matches certain known protein sequence segments by chance, these chance hits are likely to scatter in any of the three frames
FRAME SCORE = $variance_{i \in {0, 1, 2}}{S_{i}} = \frac{\sum_{i=0}^{2}(S_{i} - \bar{S})^2}{2}$

http://cpc.cbi.pku.edu.cn

### Data mining: differential expression and clustering

Differentially expressed genes
Co-expressed genes

Data mining: differentially expression calling
- identigy the genes with biological-significant difference in expression levels across samples
- differences in expression values can result from many non-biological sources(e.g. experiment error/bias)
  - the 'real' differences are the differences that can NOT be explained by the various errors introduced during the expreimental phase

1. Random errors arise from random fluctuations in the measurements
2. It could be reduced by repeating experiment many times (and get a mean value)
3. Ranodm errors could be modeled statistically by variance.

Statistical calling
1. Select a statistic whick takes the variance into account, and will rank the genes in order of supporting strength for 'differential expression'.
2. Derive the p-value for each gene, based on the NULL distribution of the statistic.
3. Choose a critical-value for the gene with p-value less than which being called as 'being statistically significant'.
$\frac{signal}{noise} = \frac{difference\, between\, group\, means}{variability\, of\, groups} = \frac{\bar{X}_{T} - \bar{X}_{C}}{SE(\bar{X}_{T} - \bar{X}_{C})}$ = t-value
  - the t-test assesses thether the means of two groups are statistically different from each other
    - take the variance into account through Standard Error (SE)
  - Need to estimate the SE correctly
    - but the correct estimation depends on prior distribution (Normal) as well as the number of replicates (>10)

Hypothesis test

|Output of statistical test| Hypothesis truth?|Hypothesis truth?|
|------|------|------|
|| $H_{1}$ (active)| $H_{0}$ (inactive)|
|Reject $H_{0}$ (active)|Hit|Type I error|
|Accept $H_{0}$ (inactive)|Type II error| Correct rejection|

Type I Error (False Positive): rejecting the null hypothesis when it is true
Type II Error (False Negative): accepting the null hypothesis when it is false

Multiple testing issue
- if more than one test is made, then the collective FP value is greater than in the single-test
  - that is, overall Type I error increases
- e.g: you checked your RNA-Seq data and found 20 significantly different genes with a 0.05 threshold on each gene, then what is the chance that you making at least one error in overall?
- Pr(making a mistake) = 0.05
- Pr(not making a mistake) = $1-0.05=0.95$
- Pr(not making any mistake) = $0.95^{20} = 0.358$
- Pr(making at least one mistake) = $1-0.358=0.642$

There is a 64.2% chance of making at least one mistake

Bonferroni Correction
- Most straightforward and plain
- For n hypothesis tests, only call p-values less than $\alpha$/n as "being significant".
  - Or, adjust the raw p-value as min(n*p, 1)
- For example, if we want to have an experiment wide Type I error rate of 0.05 when we comparing 30000 genes, we'd need p-values called as "being significant"

Type I (false positive) error rates

||# not rejected|# rejected|totals|
|------|------|------|------|
|# true H|U|V (False Positive)|$m_{0}$|
|# non-true H|T (False Negative)|S|$m_{1}$|
|totals|m-R|R|m|

- Family-wise Error Rate: $FWER = p(V \geqslant 1)$
- Per-family Error Rate: $PFER = E(V)$
- Per-comparison Error Rate: $PCER = E(V)/m$
- False Discovery Rate: $FDR = E(V/R)$
- False Positive Rate: $FPR = E(V/m_{0})$

q-value
- q-value is an measure of False Discovery Rate (FDR) 
  - Proposed by Storey et al. in 2002 and tuned for microarray analysis
- The q-value for a particular gene g is the expected proportion of false positives incurred when calling that gene g "significant"
- In contrast, the p-value for a particular gene g is the probalility that a randomly generated expression profile would be as or more extremely differentially expressed.
- Differentially expressed genes
- Co-expressed genes

Clustering: group cases (genes/samples) with similar expression pattern/levels (Unsupervised learning)
- Hierarchical Cluster, k-mean Cluster, Self-Organizing Maps (SOM), etc.

Distance measurement: how "similar" between two genes' profile
Euclidean distance (Absolution distance): $s(x_{1}, x_{2}) = \sqrt{\sum(x_{1k}^2 - x_{2k}^2)}$

Pearson distance (Correlation distance): $s(x_{1}, x_{2}) = \frac{\sum_{k=1}^{K}(x_{1k} - \bar{x_{1}})(x_{2k}-\bar{x_{2})}}{\sqrt{\sum_{k=1}^{K}(x_{1k}-\bar{x_{1}})^2 \sum_{k=1}^{K}(x_{2k}-\bar{x_{2}})^2}}$

Readings

Required

Genome Biol. 2013 Sep 10;14(9):R95. Comprehensive evaluation of differential gene expression analysis methods for RNA-seq data. Rapaport F, Khanin R, Liang Y, Pirun M, Krek A, Zumbo P, Mason CE, Socci ND, Betel D. PMID: 24020486

week 10

## Ontology and Identification of Molecular Pathways

### Ontology and Gene Ontology







OBO format
XML format
RDF-XML format

http://www.geneontology.org/

### KEGG Pathway Database

http://www.kegg.jp/kegg/

Flat format
KEGG Markup Language (KGML) format


### Pathway Identification

http://kobas.cbi.pku.edu.cn/


### An Application: common molecular pathways underlying addiction

J. Widom. Introduction to Databases. https://www.coursera.org/course/db

Readings Required

Nat Genet. 2000 May;25(1):25-9. Gene ontology: tool for the unification of biology. The Gene Ontology Consortium. PMID: 10802651

Bioinformatics. 2005 Oct 1;21(19):3787-93. Automated genome annotation and pathway identification using the KEGG Orthology (KO) as a controlled vocabulary. Mao X, Cai T, Olyarchuk JG, Wei L. PMID: 15817693

Nucleic Acids Res. 2011 Jul;39(Web Server issue):W316-22. KOBAS 2.0: a web server for annotation and identification of enriched pathways and diseases. Xie C, Mao X, Huang J, Ding Y, Wu J, Dong S, Kong L, Gao G, Li CY, Wei L. PMID: 21715386

PLoS Comput Biol. 2008 Jan;4(1):e2. Genes and (common) pathways underlying drug addiction. Li CY, Mao X, Wei L. PMID: 18179280

week 11

## 

### 

###

### European Bioinformatics Institute


Ensembl (http://ensembl.org)

IntAct (http://www.ebi.ac.uk/intact)

Clustal Omega (http://www.clustal.org/omega/)

ClustalW2 - Phylogeny (http://www.ebi.ac.uk/Tools/phylogeny/clustalw2_phylogeny/)

InterPro & InterProScan (http://www.ebi.ac.uk/interpro/)

PRIDE (http://www.ebi.ac.uk/pride/)

### UCSC Genome Browser

UCSC Genome Bioinformatics/ Browser (http://genome.ucsc.edu/)

http://genome.ucsc.edu/cgi-bin/hgTracks

ENCODE Project (http://www.genome.gov/10005107)

ENCODE data portal at UCSC (http://genome.ucsc.edu/encode/)

Neandertal Genome (http://genome.ucsc.edu/Neandertal/)

UCSC BLAT (http://genome.ucsc.edu/cgi-bin/hgBlat)

UCSC In-Silico PCR (http://genome.ucsc.edu/cgi-bin/hgPcr)

### Individual Resources

PDB (http://www.rcsb.org/)

SWISS-MODEL (http://swissmodel.expasy.org/workspace/index.php?func=modelling_simple1&userid=USERID&token=TOKEN)

SWISS-MODEL Repository (http://swissmodel.expasy.org/repository/?pid=smr01&zid=async)

I-TASSER (http://zhanglab.ccmb.med.umich.edu/I-TASSER/)

QUARK (http://zhanglab.ccmb.med.umich.edu/QUARK/)

Protein Model Portal (http://www.proteinmodelportal.org/)
meta search method which contains 6 template-based predict methods: ModWeb, M4T, SWISS-MODEL, I-TASSER, HHpred, Phyre2

CASP (http://www.predictioncenter.org/)

modENCODE (http://www.modencode.org/)

Rfam (http://rfam.sanger.ac.uk, new site http://rfam.xfam.org/)

PlantTFDB: Plant Transcription Factor Database (http://planttfdb.cbi.pku.edu.cn/)

SOAPdenovo (http://soap.genomics.org.cn/soapdenovo.html)

CNVnator (http://sv.gersteinlab.org/cnvnator/)

Three tools:
- Bioconductor (http://bioconductor.org/)
- BioPerl (http://bioperl.org/)
- BioPython (http://biopython.org)

### CBI Resources Review

http://www.cbi.pku.edu.cn/research/

week 12

# New Gene Evolution Detected by Genomic Computation

## Basic Concepts and Examples

### Case Study I: Origination of New Genes

## A Driver for Human Brain Evolution


Readings

Required

Nat Rev Genet. 2013 Sep;14(9):645-60. New genes as drivers of phenotypic evolution. Chen S, Krinsky BH, Long M. PMID: 23949544  

PLoS Biol. 2010 Oct 5;8(10). Chromosomal redistribution of male-biased genes in mammalian evolution with two bursts of gene gain on the X chromosome. Zhang YE, Vibranovski MD, Landback P, Marais GA, Long M. PMID: 20957185  

PLoS Comput Biol. 2010 Mar 26;6(3):e1000734. A human-specific de novo protein-coding gene associated with human brain functions. Li CY, Zhang Y, Wang Z, Zhang Y, Cao C, Zhang PW, Lu SJ, Li XM, Yu Q, Zheng X, Du Q, Uhl GR, Liu QR, Wei L. PMID: 20376170  

PLoS Genet. 2012 Sep;8(9):e1002942. Hominoid-specific de novo protein-coding genes originating from long non-coding RNAs. Xie C, Zhang YE, Chen JY, Liu CJ, Zhou WZ, Li Y, Zhang M, Zhang R, Wei L, Li CY. PMID: 23028352  

week 13

# Evolution Function Analysis of DNA Methyltransferase

## From Dry to Wet, an Evolutionary Story




week 14




