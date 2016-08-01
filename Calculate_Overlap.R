#	Calculate Overlap: Program that calculates the overlap count between each of the
#		TFBS ChIP-seq experiments and the 1000 randomized LD block files. This program
#		can also be modified and used to calculate the overlap count of the monocyte
#		LD blocks and the TFBS ChIP-seq experiments.
#
#	Author: Choudhury, Mudra
#
#	Purpose: To calculate the overlap between each of the randomized LD block files
#		(or monocyte LD block file) and each TFBS ChIP-seq experiment.
#
#	Usage: This program uses the ChIP-seq experiments (which must be .bed files) 
#		and the randomized LD block files titled "LDrandomization_[number].bed.txt"
#		which must be obtained by using the program called "Randomization.R".
#		It calculates the overlap between each TFBS experiment and each of the 1000
#		randomized LD block files. The overlap counts for each TFBS experiment are 
#		recorded in a file that is saved as "randomization.tfbs.overlap.bed.txt".
#		If this program is used to calculate the actual overlap between the TFBS 
#		experiments and the monocyte LD blocks, then the loop from 1:1000 can be 
#		removed and only the "Monocyte_LD_Blocks.bed.txt" file can be used.
#
# 
# This file is part of the [insert name of page] software package. [insert name of page] is free
# software: you can redistribute it and/or modify it under the terms of
# Apache Software License version 2.0 (and incorporated into this package
# distribution in the file LICENSE).
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the Apache Software License version 2.0 for 
# details.#
# Copyright Mudra Choudhury and Stephen Ramsey, Oregon State University
# 2016.02
#


library(rtracklayer)
library(stringr)
library(GenomicRanges)

#save all file names
myfilenames <- dir()



# Function which takes each LD block randomization file named "LDrandomization[number].bed.txt
#	which was created in the "Randomization.R" program and calls a function to calculate the
#	overlap between the randomized LD block file and each TFBS experiment.
ld.randomizer <- function(){
    for (n in 1:1000){
    print(n)
	random.snp.filename <- sprintf("LDrandomization_%d.bed.txt",n)
	randomization <- multi.count.peaks(random.snp.filename, n) 
    }
}




# Function which takes in the randomized LD block file, calls the function "count.peaks" in order
#	to calculate the overlap between each ChIP-seq experiment and the randomized LD block file, 
#	then saves the overlap counts for each ChIP-seq experiment and all 1000 randomized LD block 
#	files in a file called "randomization.tfbs.overlap.bed.txt".
multi.count.peaks <- function(snps.file, n){

   tf.bed.files <- grep(".bed$",myfilenames)
   tf.overlapped.files <- grep("^overlap.bed", myfilenames)
   all.tf.files <- setdiff(tf.bed.files, tf.overlapped.files)
   all.tf.file.names <- myfilenames[all.tf.files]
   retvec <- sapply(all.tf.file.names, count.peaks, snps.file, n)
   retvec.data <- as.data.frame(retvec)
   retvec.data1 <- str_split_fixed(retvec.data$retvec, "\t", 1)
   retvec.final <- cbind(all.tf.file.names, retvec.data1)
   
   if (n==1){
   tfbs.overlap.filename <- sprintf("randomization.tfbs.overlap.bed.txt")
   write.table(retvec.final, file=tfbs.overlap.filename, quote= FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
   }
   
   else if (n==2){
   	old.tfbs.overlap.file <- read.table("randomization.tfbs.overlap.bed.txt")
   	new.tfbs.overlap.file <- cbind(all.tf.file.names, old.tfbs.overlap.file, retvec.data1)
   write.table(new.tfbs.overlap.file, file="randomization.tfbs.overlap.bed.txt", quote= FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
   }
   
   else{
   		old.tfbs.overlap.file <- read.table("randomization.tfbs.overlap.bed.txt")
   	new.tfbs.overlap.file <- cbind(old.tfbs.overlap.file, retvec.data1)
   write.table(new.tfbs.overlap.file, file="randomization.tfbs.overlap.bed.txt", quote= FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
   }
   
}



# Function that takes in a TFBS file and a randomized LD block file as an argument and calculates
#	the overlap between the two files. It returns the overlap count and can also be used to obtain
#	a file of the overlap locations.


count.peaks <- function(tfbs.file, snps.file, n){ 
	tfbs.table <- read.table(tfbs.file) #Read in the tfbs file
	tfbs <- tfbs.table[,c("V1", "V2", "V3")] #only keep chromosome, and both boundaries
	write.table(tfbs, file=tfbs.file, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE) #save this where original file was (because other columns are irrelevant)
	tfbs <- import.bed(tfbs.file) #now you are able to import the file as a bed file
	tfbs <- unique(tfbs) 
	snps.original <- read.table(snps.file) #import ld blocks
	colnames(snps.original) <- c("chr", "start", "end")
	snps <- makeGRangesFromDataFrame(snps.original)

	#convert both data into ranges lists

	tfbs.ranges <- as(tfbs, "RangesList")
	tfbs.chroms <- names(tfbs.ranges)
	snps.ranges <- as(snps, "RangesList")
	snps.chroms <- names(snps.ranges)

	remove2 <- c("chrX", "chrY")
	snps.chroms <- snps.chroms[! snps.chroms %in% remove2]
	tfbs.chroms <- tfbs.chroms[! tfbs.chroms %in% remove2]

	#match up chromosome numbers

	inds.match.snps <- match(snps.chroms,tfbs.chroms)
	inds.match.tfbs <- match(tfbs.chroms, snps.chroms)
	inds.match.tfbs <- inds.match.tfbs[! is.na(inds.match.tfbs)]

	#make sure both have the same number of chromosomes being compared

	snps.chroms[inds.match.tfbs]
	tfbs.chroms[inds.match.snps]

	chroms.unique <- unique(snps.chroms[inds.match.tfbs])

	#make sure the order in which the chromosomes are organized are the same.

	inds.match.snps <- match(chroms.unique, snps.chroms)
	inds.match.tfbs <- match(chroms.unique, tfbs.chroms)

	remove <- c(23, 24)
	inds.match.snps <- inds.match.snps[! inds.match.snps %in% remove]
	inds.match.tfbs <- inds.match.tfbs[! inds.match.tfbs %in% remove]


	snps.chroms[inds.match.snps]
	tfbs.chroms[inds.match.tfbs]
	tfbs.ranges[inds.match.tfbs]

	names(tfbs.ranges[inds.match.tfbs])
	names(snps.ranges[inds.match.snps])

	snps1 <- snps.ranges[inds.match.snps]
	tfbs1 <- tfbs.ranges[inds.match.tfbs]

	#find overlapping regions

	overlap <- findOverlaps(tfbs1, snps1)
	
	overlapped.ranges <- ranges(overlap, tfbs1, snps1)

	overlap.file <- as.data.frame(overlapped.ranges)


	#count the total number of overlapping regions.

	count <- length(queryHits(overlap))

	#save overlapping regions as a .bed file if you would like to in next command.

	#write.table(overlap.file, file=filename, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

return(count)

}
