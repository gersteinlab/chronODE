#!/usr/bin/env Rscript


# A script for classifying curved fitted by chronODE 
# as accelerators, decelerators, or switchers

#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c ("-p", "--paramsfile"),
               help = "Parameters file to use as input" ),

  make_option( c ("-f", "--valuesfile"),
               help = "Values file to use as input" ),

  make_option( c ("-o", "--outfile"),
               help = "Name of output file." ),

  make_option( c ("-s", "--start"), type = 'numeric',
               help = "Start time (in days)" ),

  make_option( c ("-e", "--end"), type = 'numeric',
               help = "End time (in days)" )
)

parser <- OptionParser(
  usage = "%prog [options] file", 
  option_list=option_list
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#********
# BEGIN *
#********

# # 0. debugging options
# setwd("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/filtering.classification/mse.classification/work/b8/3abe11bb40598b5989cd3a9f1c9f48")
# opt <- list(paramsfile = "mse.param.tsv", valuesfile = "mseClass.fitted.tsv",
#             start = 10.5, end = 21)

# 1. read chronode parameters
param <- read.csv(opt$paramsfile, sep="\t")


# 2. compute FDR for MSE permutation p-value
param$fdr <- p.adjust(param$p_value, method = "BH")


# 3. read fitted values 
values <- read.csv(opt$valuesfile, sep="\t")


# 4. make sure param and values have same order of elements
stopifnot(identical(param[, 1], values[, 1]))


# 5. prepare vectors to store
# the switching, saturation and minimum times
switching_time <- c()
saturation_time <- c()
minimum_time <- c()


# 6. Compute, for every element, the switching, minimum, and saturation times
for (i in 1:nrow(param)) {
    
  b = param[i, "b"]
  k = param[i, "k"]
  
  y0 <- values[i, 2]
  c = log(y0) - log(1 - (y0/b)) - (k * opt$start)
  C = exp(c)
  
  switching_time <- c(switching_time, (log(b/C)/k))
  saturation_time <- c(saturation_time, (log((0.99*(b^2))/((b*C)-(0.99*b*C)))/k))
  minimum_time <- c(minimum_time, (log(((10^-16)*(b^2))/((b*C)-((10^-16)*b*C)))/k))
  
}

param$switching_time <- switching_time
param$saturation_time <- saturation_time
param$minimum_time <- minimum_time


# 7. Compute new label for each element
# based on the switching time
label <- c()

for (i in 1:nrow(param)) {
  
  if (param[i, "switching_time"] > opt$end) {
    
    label <- c(label, "accelerator")
    
  } else if (param[i, "switching_time"] < opt$start) {
      
      label <- c(label, "decelerator")
    
  } else {
      
      label <- c(label, "switcher")
      
    }
    
}

param$label <- label


# Write output table!
write.table(param, file = opt$outfile, sep = '\t', quote = F, row.names = F, col.names = T)





