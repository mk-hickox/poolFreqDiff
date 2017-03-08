# poolFreqDiff

Scripts for the analysis of allele frequency differences in genomic data from stratified populations

These scripts should be useful in the analysis of e.g. pool-seq data from experimental evolution lines.

All scripts have been tested for python 2.7 and R 3.3.2 in Ubuntu 14.

One step of this analysis requires that the programming language [R](https://www.r-project.org/) be installed on the computer

The directory includes the following scripts:


## 1) poolFreqDiffTest_QBGLM.py

This script takes a PoPoolation/PoPoolation2-style ".sync" file and outputs an R script that can be run as a batch job.

For a help message and some details type

	poolFreqDiffTest_QBGLM.py -h


## 2) poolFreqDiffTest.R

This is a collection of R functions that are needed to run the R script (see below).
This file should be kept in the same directory as poolFreqDiffTest_QBGLM.py


## 3) G_test.R

This is a function written in R that performs the G-test as described in: Sokal and Rohlf 1969,1981 "Biometry", W.H. Freeman and Company, San Francisco
This file should be kept in the same directory as poolFreqDiffTest_QBGLM.py


## 4) Example

The directory also includes a short test file which is a standard .sync file as implemented in: 

PoPoolation (Kofler et al., 2011 PLoS ONE 6:e15925) and PoPoolation2 (Kofler et al., 2011 Bioinformatics 27:3434 - 3436)

There are two SNPs in this file (line 2 [2:92271] and line 7 [2:92276]), one of which (line 7 [2:92276]) has a highly significant allele frequency difference between two treatments

To run the analysis, within the directory of scripts simply type:

	python ./poolFreqDiffTest_QBGLM.py -filename test.sync -npops 4 -nlevels 2 -n 40 -mincnt 2 -minc 10 -maxc 100 -rescale nr -zeroes 1 > test.rin

This will create an Rscript "test.rin" which can then be run:

	Rscript test.rin > test.rout

This will create an output table which includes all the same information as in the original .sync file and an additional column with the treatment effect p-value.

