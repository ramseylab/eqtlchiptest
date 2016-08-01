#	Randomize LD Blocks: Program that creates 1000 files of randomized LD blocks 
#		according to the background model file provided. This script uses the original LD 
#		blocks and shuffles them to nonoverlapping locations in the background model
#		while preserving their widths. The 1000 randomized files can be used in the
#		script "Calculate_Overlap.R" to calculate the overlap between each randomized LD
#		block and each TFBS ChIP-seq experiment.
#
#	Author: Choudhury, Mudra
#
#	Purpose: To create 1000 files of randomized LD blocks according to a background 
#		model. These will later used for statistical analysis.
#
#	Usage: The file of original LD blocks "Monocyte_LD_Blocks.bed.txt" will be needed. The 
#		background model file "Background_Model.txt" file is needed to restrict the
#		randomization to the background model. The output will be 1000 randomized LD block
#		files called "LDrandomization_[number].bed.txt".
#
# 
# This file is part of the [insert page name] software package. [Insert page name] is free
# software: you can redistribute it and/or modify it under the terms of
# Apache Software License version 2.0 (and incorporated into this package
# distribution in the file LICENSE).
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the Apache Software License version 2.0 for 
# details.
#
# Copyright Mudra Choudhury and Stephen Ramsey, Oregon State University
# 2016.02
#

library(rtracklayer)
library(stringr)
library(GenomicRanges)


#This function takes in the original monocyte LD blocks "LD_Blocks.bed.txt" and shuffles the LD blocks to create a randomized files.
ld.randomizer <- function(snps.file){
	x <- import.bed(snps.file, genome="hg19") #Load in "LD_Blocks.bed.txt"
	for (n in 1:1000){
		genome_file <- tempfile()
		
		#save widths of original LD blocks
		writeLines(paste(names(seqlengths(x)), seqlengths(x), sep="\t"), genome_file) 
		bed_file <- tempfile()
    	x <- x[order(width(x), decreasing=TRUE)]
    	writeLines(paste(as.character(seqnames(x)), start(x), end(x), sep="\t"), bed_file)
    	shuffle_file <- tempfile()
    	
    	#shuffle the LD blocks to nonoverlapping locations according to the background model
    	system(paste("shuffleBed -maxTries 10000 -noOverlapping -incl Background_Model.txt -i", bed_file, "-g", genome_file, ">", shuffle_file))
    	shuffled <- import.bed(shuffle_file)
    	seqlevels(shuffled, force=TRUE) <- seqlevels(x)
    	seqlengths(shuffled) <- seqlengths(x)
		unlink(c(genome_file, bed_file, shuffle_file))
		snp.data.shuffled <- as.data.frame(shuffled)
		
		#Shift LD blocks to the right by half of their width
		half.width <- (snp.data.shuffled$width/2)
		snp.data.shuffled$start <- (snp.data.shuffled$start)-half.width
		snp.data.shuffled$end <- (snp.data.shuffled$start + snp.data.shuffled$width)
		shuffled2 <- makeGRangesFromDataFrame(snp.data.shuffled)
		shuffled <- trim(shuffled2)
		random.snp.file <- shuffled
		
		#save the file of randomized LD blocks
		random.snp.filename <- sprintf("LDrandomization_%d.bed.txt",n)
		export.bed(random.snp.file, random.snp.filename)
		
	}
	
}