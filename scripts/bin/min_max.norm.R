#!/usr/bin/env Rscript


#*************
# author: BB |
#*************


#************
# LIBRARIES *
#************

library(data.table)
library(optparse)


#*****************
# OPTION PARSING *
#*****************

option_list <- list(
  
  make_option( c ( "-i", "--input_matrix" ),
               help = "The matrix obtained after normalization.")
)

parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options



#************
# FUNCTIONS *
#************

# generalized form of the min-max normalization
# where a is the desired new min value (e.g. 0)
# and b is the desired new max value (e.g. 1)
normalize <- function(x, a, b, na.rm = TRUE) {
  return( (((b - a)*(x - min(x))) / (max(x)-min(x))) + a )
}


#********
# BEGIN *
#********


# # 0. debugging options
# setwd("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/work/f4/cfb70f2e7b0559b78cdd4928143dd4")
# opt <- list()
# opt$input_matrix <- "forebrain.avg.2.tsv"

# 1. read input matrix
m <- fread(opt$input_matrix, data.table = F, h=T, stringsAsFactors = F)


# 2. rescale global shift matrix to min-max 0-1 range
m.minmax.0_1 <- as.data.frame(t(apply(as.data.frame(m[, 2:ncol(m)]), 1, function(x){x <- normalize(x = x, a = 0.00001, b = 1)})))
m.minmax.0_1[, ncol(m.minmax.0_1)+1] <- m[, 1]
colnames(m.minmax.0_1)[ncol(m.minmax.0_1)] <- colnames(m)[1]
m.minmax.0_1 <- m.minmax.0_1[, colnames(m)]


# 3. rescale global shift matrix to min-max 1-2 range
m.minmax.1_2 <- as.data.frame(t(apply(as.data.frame(m[, 2:ncol(m)]), 1, function(x){x <- normalize(x = x, a = 1, b = 2)})))
m.minmax.1_2[, ncol(m.minmax.1_2)+1] <- m[, 1]
colnames(m.minmax.1_2)[ncol(m.minmax.1_2)] <- colnames(m)[1]
m.minmax.1_2 <- m.minmax.1_2[, colnames(m)]


# 4. save outputs
write.table(m.minmax.0_1, file = "minmax.0_1.tsv", row.names = F, col.names = T, sep = "\t", quote = F)
write.table(m.minmax.1_2, file = "minmax.1_2.tsv", row.names = F, col.names = T, sep = "\t", quote = F)




