#--------------------------------------------------#
# This script will read a normal .sync file        #
# and plot the allele frequencies of the SNPs.     #
#                                                  #
# Author: R. Axe. W. Wiberg                        #
# Created: March 2017                              #
# Last modified: March 2017                        #
#--------------------------------------------------#

# N.B. IF THERE ARE MORE THAN 50 SNPs THIS SCRIPT WILL 
# PLOT A RANDOM SET OF 50.

# N.B. THIS SCRIPT ASSUMES:


# Load libraries
library("optparse")
library("methods")
# Set maximum number of decimal positions to a
# high number.
options(scipen = 999)

option_list = list(
  make_option(c("-f","--file"), type="numeric", default=0, 
              help="input .sync file", 
              metavar="numeric")
  )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print_help(opt_parser)

read.table(opt$file, sep="\t",header = FALSE)