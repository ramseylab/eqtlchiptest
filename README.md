eqtlchiptest: Software for statistical testing of the overlap between expression
quantitative trait loci (eQTLs) and genome-wide transcription factor location
data.  This software accompanies the manuscript "Identifying cell type-specific
transcription factors by integrating ChIP-seq and eQTL data--application to
monocyte gene regulation", by Mudra Choudhury and Stephen Ramsey, which has been
submitted to the Journal of Bioinformatics and Computational Biology.

Author:  Mudra Choudhury, Oregon State University

From the lab of:  Stephen Ramsey, Oregon State University (lab.saramsey.org)

Date:  August 1, 2016

This software is distributed under the Apache Software License 2.0.
Please see the file LICENSE for details on the software licensing
agreement.

Usage notes: There are two R scripts in the subdirectory "R", that comprise this
software release. The scripts load tab-delimited text datafiles, examples of
which are given in the subdirectory "data".  Each of the example datafiles
contains the first ten lines of the actual datafile used in the analysis. The
example datafiles have Unix line termination. The two R scripts are used in
this order:

(1) Randomization.R
(2) Calculate_Overlap.R

The "Randomization.R" script generates a file (UCSC BED format) of eQTLs that
are randomly placed in the genome, within a set of allowed regions defined by
the file "Background_Model.txt" (UCSC BED format). The sizes of the randomly
placed eQTL regions are defined by the input file "Monocyte_LD_Blocks.bed".

The "Calculate_Overlap.R" script compares an input file of eQTL regions (in UCSC
BED format) with an input file of ChIP-seq peak regions (in UCSC BED format) to
compute the number of ChIP-seq peaks whose centers are within any of the eQTL
regions. The sex chromosomes X and Y are excluded from the analysis (this
exclusion is because when the authors used the SNAP tool to map eSNPs to proxy
SNPs in order to define the eQTLs, proxy SNPs on the sex chromosomes were not
returned by SNAP (see the above-referenced manuscript for more details).



