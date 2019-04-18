---
layout: page
permalink: /tools/
description: "Other tools for analysis high-throughput experiments."
modified: 2013-09-11
tags: [Bowtie, Tophat, Cufflinks, CummeRbund]
---

# Bowtie: ultrafast short read alignment

[Bowtie](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters, and particularly good at aligning to relatively long (e.g. mammalian) genomes. Bowtie 2 indexes the genome with an FM Index to keep its memory footprint small: for the human genome, its memory footprint is typically around 3.2 GB. Bowtie 2 supports gapped, local, and paired-end alignment modes.

# TopHat: alignment of short RNA-Seq reads

[TopHat](http://ccb.jhu.edu/software/tophat/index.shtml) is a fast splice junction mapper for RNA-Seq reads. It aligns RNA-Seq reads to mammalian-sized genomes using the ultra high-throughput short read aligner Bowtie, and then analyzes the mapping results to identify splice junctions between exons. 

TopHat is a collaborative effort among Daehwan Kim and Steven Salzberg in the Center for Computational Biology at Johns Hopkins University, and Cole Trapnell in the Genome Sciences Department at the University of Washington. TopHat was originally developed by Cole Trapnell at the Center for Bioinformatics and Computational Biology at the University of Maryland, College Park.

TopHat is provided under the OSI-approved Artistic License 2.0.

# Cufflinks: transcript assembly and differential expression for RNA-Seq

[Cufflinks](http://cufflinks.cbcb.umd.edu/) assembles transcripts, estimates their abundances, and tests for differential expression and regulation in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundances of these transcripts based on how many reads support each one, taking into account biases in library preparation protocols. 

Cufflinks was originally written by Cole Trapnell while working with Lior Pachter at UC Berkeley. The project was a collaborative effort between the Pachter lab, Steven Salzberg's computational genomics group at the Institute of Genetic Medicine at Johns Hopkins University, and Barbara Wold's lab at Caltech. 

Cufflinks is provided under the OSI-approved Boost License.

# CummeRbund: visualization of RNA-Seq differential analysis

[CummeRbund](http://compbio.mit.edu/cummeRbund/) is an R package that is designed to aid and simplify the task of analyzing Cufflinks RNA-Seq output.
CummeRbund was written by Loyal Goff while in Manolis Kellis's group at MIT's Computer Science and Artificial Intelligence Laboratory, and John Rinn's Lab at the Harvard University department of Stem Cells and Regenerative Medicine.

CummeRbund is provided under the OSI-approved Artistic License 2.0.
