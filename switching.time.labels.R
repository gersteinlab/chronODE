
#!/usr/bin/env Rscript

.libPaths("/gpfs/gibbs/pi/gerstein/bb926/R/libraries/")


#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c ("-r", "--region"),
               help = "The brain subregion (fore/mid/hindbrain)" ),

  make_option( c ("-g", "--group"),
               help = "Either upreg or downreg." ),
  
  make_option( c ("-o", "--outFile"),
               help = "Name of output file." )

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

# 0. debugging options
# opt <- list()
# opt$region <- "forebrain"
# opt$group <- "upreg"

# 1. read dataframe with set of 
# cCREs classified based on concave/convex approach by Ke
m <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/open.chromatin/clustering/convex.concave/", 
                     opt$region, "_label.csv"),
              sep="\t")

# 2. filter cCREs based on the group
if (opt$group == "upreg") {
  
  m <- m[m$group != "downreg", ]
  
} else {
  
  m <- m[m$group == "downreg", ]
  
}

# 3. read actual values reconstructed
# by the ODE pipeline
if (opt$group == "upreg") {
  
  values.lu <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/open.chromatin/values/", 
                               opt$region, ".late_upreg.jan5.tsv"),
                        sep="\t")
  values.eu <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/open.chromatin/values/", 
                               opt$region, ".early_upreg.jan5.tsv"),
                        sep="\t")
  values <- rbind(values.lu, values.eu)
  
} else {
  
  values <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse/ODE/open.chromatin/values/", 
                               opt$region, ".downreg.jan5.tsv"),
                        sep="\t")
}

# 4. prepare vectors to store
# values of Ce^(kt) at 10.5 and 21 days,
# and the switching and saturation times
cekt_10 <- c()
cekt_21 <- c()
switching_time <- c()
saturation_time <- c()
minimum_time <- c()

# 5. Compute, for every cCRE, Ce^(kt) at 10.5 and 21 days
# and the switching and saturation times

for (i in 1:nrow(m)) {
  
  cCRE_id = m$cCRE_id[i]
  
  b = m[m$cCRE_id == cCRE_id, "b"]
  k= m[m$cCRE_id == cCRE_id, "k"]
  
  y0 <- values[values$cCRE_id == cCRE_id, 2]
  c = log(abs(y0)) - log(abs(1 - (y0/b))) - k*10.5
  C = exp(c)
  
  cekt_10 <- c(cekt_10, C*exp(k*10.5))
  cekt_21 <- c(cekt_21, C*exp(k*21))
  switching_time <- c(switching_time, (log(b/C)/k))
  saturation_time <- c(saturation_time, (log((0.99*(b^2))/((b*C)-(0.99*b*C)))/k))
  minimum_time <- c(minimum_time, (log(((10^-16)*(b^2))/((b*C)-((10^-16)*b*C)))/k))
  
}

sub.m <- m
sub.m$cekt_10 <- cekt_10
sub.m$cekt_21 <- cekt_21
sub.m$ratio_10 <- sub.m$cekt_10 / sub.m$b
sub.m$ratio_21 <- sub.m$cekt_21 / sub.m$b
sub.m$switching_time <- switching_time
sub.m$saturation_time <- saturation_time
sub.m$minimum_time <- minimum_time

# 6. keep only cCREs that have consistent 
# increase (upreg) or decrease (downreg) 
# in the Ce^(kt) / b ratio between 10.5 and 21 days
if (opt$group == "upreg") {
  
  sub.m <- sub.m[sub.m$ratio_10 < sub.m$ratio_21, ]
  
} else {
  
  sub.m <- sub.m[sub.m$ratio_10 > sub.m$ratio_21, ]
  
}

# 7. Compute new label for each cCRE
# based on the switching time
label_switching_time <- c()

for (i in 1:nrow(sub.m)) {
  
  if (sub.m[i, "switching_time"] > 21) {
    
    label_switching_time <- c(label_switching_time, "accelerator")
    
  } else if (sub.m[i, "switching_time"] < 10.5) {
      
      label_switching_time <- c(label_switching_time, "decelerator")
    
  } else {
      
      label_switching_time <- c(label_switching_time, "switcher")
      
    }
    
}

sub.m$label_switching_time <- label_switching_time
colnames(sub.m)[colnames(sub.m) == "label"] <- "label_convex_concave"

table(sub.m[, c("label_convex_concave", "label_switching_time")])

# 8. save output matrix
write.table(sub.m, file = opt$outFile,
            sep = '\t', quote = F, row.names = F, col.names = T)


# # plotting examples

# test <- values[values$cCRE_id == "EM10D1042940", 2:ncol(values)]
# t.test <- t(test)
# colnames(t.test) <- "value"
# t.test <- as.data.frame(t.test)
# t.test$tp <- gsub("_oc", "", rownames(t.test))
# t.test$tp <- gsub("X", "", t.test$tp)
# t.test$tp <- as.numeric(t.test$tp)
# 
# ggplot(t.test, aes(x=tp, y=value)) +geom_line()


