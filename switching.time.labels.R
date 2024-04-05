

#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c ("-p", "--paramsfile"),
               help = "Parameters file to use as input" ),

  make_option( c ("-v", "--valuesfile"),
               help = "Values file to use as input" ),

  make_option( c ("-o", "--outfile"),
               help = "Name of output file." ),

  make_option( c ("-s", "--start"), default = 10.5,
               help = "Start time (in days)" ),

  make_option( c ("-e", "--end"), default = 21,
               help = "End time (in days)" ),

  make_option( c ("-g", "--group"), default = "none",
               help = "Pattern filter (increasing or decreasing)" )
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

# 1. read dataframe with set of 
# cCREs classified based on concave/convex approach by Ke
m <- read.csv(opt$paramsfile, sep="\t")

# 3. read actual values reconstructed
# by the ODE pipeline
values <- read.csv(opt$valuesfile, sep="\t")

# 4. prepare vectors to store
# values of Ce^(kt) at start and end times,
# and the switching and saturation times
cekt_start <- c()
cekt_end <- c()
switching_time <- c()
saturation_time <- c()
minimum_time <- c()

# 5. Compute, for every cCRE, Ce^(kt) at start and end times
# and the switching and saturation times

for (i in 1:nrow(m)) {
  
  cCRE_id = m$cCRE_id[i]
  
  b = m[m$cCRE_id == cCRE_id, "b"]
  k= m[m$cCRE_id == cCRE_id, "k"]
  
  y0 <- values[values$cCRE_id == cCRE_id, 2]
  c = log(abs(y0)) - log(abs(1 - (y0/b))) - k * opt$start
  C = exp(c)
  
  cekt_start <- c(cekt_start, C*exp(k*opt$start))
  cekt_end <- c(cekt_end, C*exp(k*opt$end))
  switching_time <- c(switching_time, (log(b/C)/k))
  saturation_time <- c(saturation_time, (log((0.99*(b^2))/((b*C)-(0.99*b*C)))/k))
  minimum_time <- c(minimum_time, (log(((10^-16)*(b^2))/((b*C)-((10^-16)*b*C)))/k))
  
}

sub.m <- m
sub.m$cekt_start <- cekt_start
sub.m$cekt_end <- cekt_end
sub.m$ratio_start <- sub.m$cekt_start / sub.m$b
sub.m$ratio_end <- sub.m$cekt_end / sub.m$b
sub.m$switching_time <- switching_time
sub.m$saturation_time <- saturation_time
sub.m$minimum_time <- minimum_time

# 6. optionally keep only cCREs that have consistent 
# increasing or decreasing
# in the Ce^(kt) / b ratio between start and end times
if (opt$group == "increasing") {
  
  sub.m <- sub.m[sub.m$ratio_start < sub.m$ratio_end, ]
  
} else if (opt$group == "decreasing") {
  
  sub.m <- sub.m[sub.m$ratio_start > sub.m$ratio_end, ]
  
}

# 7. Compute new label for each cCRE
# based on the switching time
label <- c()

for (i in 1:nrow(sub.m)) {
  
  if (sub.m[i, "switching_time"] > opt$end) {
    
    label <- c(label, "accelerator")
    
  } else if (sub.m[i, "switching_time"] < opt$start) {
      
      label <- c(label, "decelerator")
    
  } else {
      
      label <- c(label, "switcher")
      
    }
    
}

sub.m$label <- label

# 8. save output matrix
sub.m$cekt_start <- NULL
sub.m$cekt_end <- NULL
sub.m$ratio_start <- NULL
sub.m$ratio_end <- NULL

write.table(sub.m, file = opt$outfile,
            sep = '\t', quote = F, row.names = F, col.names = T)





