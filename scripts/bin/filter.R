#!/usr/bin/env Rscript


#*************
# author: BB |
#*************


#************
# LIBRARIES *
#************

library(optparse)
library(data.table)


#*****************
# OPTION PARSING *
#*****************

option_list <- list(
  
  make_option( c ( "-o", "--orig" ), help = "The original values matrix."),
  make_option( c ( "-f", "--fitted" ), help = "The fitted matrix."),
  make_option( c ( "-r", "--rescfitted" ), help = "The fitted matrix rescaled to original range."),
  make_option( c ( "-d", "--deriv" ), help = "The derivatives matrix."),
  make_option( c ( "-p", "--par" ), help = "The parameters matrix."),
  make_option( c ( "-k", "--keep" ), help = "The list of elements to keep.")
  
)

parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#************
# FUNCTIONS *
#************

read_fun <- function(myfile, tokeep) {
  
  df <- fread(myfile, data.table = F, stringsAsFactors = F, header = T)
  rownames(df) <- df[, 1]
  df <- df[complete.cases(df), ]
  df <- df[ order(row.names(df)), ]
  
  return(df)
  
}


filter_fun <- function(df, x) {
  
  df <- df[rownames(df) %in% x, ]
  df <- df[x, ]
  
  return(df)
  
}



#********
# BEGIN *
#********


# 0. debugging options
# opt <- list(orig = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/matrices/open.chromatin/ODE.input/forebrain/forebrain.dynamic.global.shift.tsv",
#             fitted = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/ODE.fitting/results/forebrain/oc/forebrain.oc.8_out.values.tsv",
#             rescfitted = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/ODE.fitting/results/forebrain/oc/forebrain.oc.8_out.rescaled.values.tsv",
#             deriv = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/ODE.fitting/results/forebrain/oc/forebrain.oc.8_out.derivatives.tsv",
#             par = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/ODE.fitting/results/forebrain/oc/forebrain.oc.8_out.parameters.tsv",
#             keep = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/filtering.classification/reproducible.fitting/forebrain.oc.tokeep.tsv")

# 1. read list of elements to keep (i.e., those that show reproducible fitting)
tokeep <- fread(opt$keep, data.table = F, header = T, stringsAsFactors = F)
tokeep <- tokeep[, 1]

# 2. read dataframes & remove NAs
r_orig <- read_fun(opt$orig)
r_fitted <- read_fun(opt$fitted)
r_rescfitted <- read_fun(opt$rescfitted)
r_deriv <- read_fun(opt$deriv)
r_par <- read_fun(opt$par)

# 3. ensure that fitted, rescfitted, and deriv have same set of entries
stopifnot(identical(rownames(r_fitted), rownames(r_rescfitted)))
stopifnot(identical(rownames(r_fitted), rownames(r_deriv)))

# 4. final list of elements needs to be an intersection 
# of tokeep and rownames of r_fitted
# this is because there are some entries with reproducible solutions that have NAs in their fitted vector
tokeep2 <- intersect(tokeep, rownames(r_fitted))

# 5. filter dataframes so that they only have the rownames in tokeep2
f_orig <- filter_fun(r_orig, x = tokeep2)
f_fitted <- filter_fun(r_fitted, x = tokeep2)
f_rescfitted <- filter_fun(r_rescfitted, x = tokeep2)
f_deriv <- filter_fun(r_deriv, x = tokeep2)
f_par <- filter_fun(r_par, x = tokeep2)

# 6. double-check row order of all dataframes
stopifnot(identical(rownames(f_orig), rownames(f_fitted)))
stopifnot(identical(rownames(f_orig), rownames(f_rescfitted)))
stopifnot(identical(rownames(f_orig), rownames(f_deriv)))
stopifnot(identical(rownames(f_orig), rownames(f_par)))

# 4. save dataframes
write.table(f_orig, file = "filtered.orig.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)

write.table(f_fitted, file = "filtered.fitted.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)

write.table(f_rescfitted, file = "filtered.rescfitted.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)

write.table(f_deriv, file = "filtered.deriv.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)

write.table(f_par, file = "filtered.par.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)




