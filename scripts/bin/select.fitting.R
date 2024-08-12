#!/usr/bin/env Rscript


#*************
# author: BB |
#*************


#************
# LIBRARIES *
#************

library(optparse)


#*****************
# OPTION PARSING *
#*****************

option_list <- list(
  
  make_option( c ( "--p_0_1" ), help = "The 0-1 parameters matrix."),
  make_option( c ( "--p_1_2" ), help = "The 1-2 parameters matrix."),
  make_option( c ( "--v_0_1" ), help = "The 0-1 values matrix."),
  make_option( c ( "--v_1_2" ), help = "The 1-2 values matrix."),
  make_option( c ( "--d_0_1" ), help = "The 0-1 derivatives matrix."),
  make_option( c ( "--d_1_2" ), help = "The 1-2 derivatives matrix."),
  make_option( c ( "--raw_values" ), help = "The raw values matrix.")
  
)

parser <- OptionParser(
  usage = "%prog [options]", 
  option_list=option_list)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options


#************
# FUNCTIONS *
#************

rescale_fun <- function(min_range, max_range, x, y) {
  
  min_val = min(x)
  max_val = max(x)
  shifted_y = (((y - min_range)*(max_val - min_val))/(max_range-min_range)) + min_val
  
  return(shifted_y)
  
}

#********
# BEGIN *
#********


# # 0. debugging options
# setwd("/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/BB.nf.test/work/83/80f31e106535e33fa31eff4e810a1c")
# opt <- list(p_0_1 = "parameters.0_1.tsv",
#             p_1_2 = "parameters.1_2.tsv",
#             v_0_1 = "values.0_1.tsv",
#             v_1_2 = "values.1_2.tsv",
#             d_0_1 = "derivatives.0_1.tsv",
#             d_1_2 = "derivatives.1_2.tsv",
#             raw_values = "forebrain.dynamic.global.shift.6.cCREs.2.tsv")


# 1. read in values
# parameters
m_0_1_parameters <- read.table(opt$p_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_parameters <- read.table(opt$p_1_2, header = T, sep="\t", stringsAsFactors = F)

# values
m_0_1_values <- read.table(opt$v_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_values <- read.table(opt$v_1_2, header = T, sep="\t", stringsAsFactors = F)

# derivatives
m_0_1_derivatives <- read.table(opt$d_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_derivatives <- read.table(opt$d_1_2, header = T, sep="\t", stringsAsFactors = F)

# raw values 
m_raw_values <- read.table(opt$raw_values, header = T, sep="\t", stringsAsFactors = F)

# compute number of time points
tp = ncol(m_0_1_values) - 1

# 2. select best fitting for each element
m_parameters <- data.frame(matrix(ncol=ncol(m_0_1_parameters),nrow=0, dimnames=list(NULL, colnames(m_0_1_parameters))), stringsAsFactors = F)
m_values <- data.frame(matrix(ncol=ncol(m_0_1_values),nrow=0, dimnames=list(NULL, colnames(m_0_1_values))), stringsAsFactors = F)
m_derivatives <- data.frame(matrix(ncol=ncol(m_0_1_derivatives),nrow=0, dimnames=list(NULL, colnames(m_0_1_derivatives))), stringsAsFactors = F)
m_rescaled_values <- data.frame(matrix(ncol=ncol(m_0_1_values),nrow=0, dimnames=list(NULL, colnames(m_0_1_values))), stringsAsFactors = F)
range <- c()
rescaled_b <- c()

for ( i in 1:nrow(m_1_2_parameters) ) {
  
  if (is.na(m_1_2_parameters[i, "MSE"]) & (is.na(m_0_1_parameters[i, "MSE"]))) {

    
    m_parameters <- rbind(m_parameters, data.frame(t(setNames(c(m_1_2_parameters[i, 1], rep(NA, 3)), names(m_parameters)))))
    m_values <- rbind(m_values, data.frame(t(setNames(c(m_1_2_values[i, 1], rep(NA, tp)), names(m_values)))))
    m_derivatives <- rbind(m_derivatives, data.frame(t(setNames(c(m_1_2_derivatives[i, 1], rep(NA, tp)), names(m_derivatives)))))
    range <- c(range, NA)
    m_rescaled_values <- rbind(m_rescaled_values, data.frame(t(setNames(c(m_1_2_values[i, 1], rep(NA, tp)), names(m_rescaled_values)))))
    rescaled_b <- c(rescaled_b, NA)
    
  } else {
    
    if (is.na(m_1_2_parameters[i, "MSE"]) & !(is.na(m_0_1_parameters[i, "MSE"]))) {
      
      sel.parameters <- m_0_1_parameters[i, ]
      sel.values <- m_0_1_values[i, ]
      sel.derivatives <- m_0_1_derivatives[i, ]
      sel.range <- "0-1"
      min_range = 1e-5
      max_range = 1
      
    } else if (is.na(m_0_1_parameters[i, "MSE"]) & !(is.na(m_1_2_parameters[i, "MSE"]))) {
      
      sel.parameters <- m_1_2_parameters[i, ]
      sel.values <- m_1_2_values[i, ]
      sel.derivatives <- m_1_2_derivatives[i, ] 
      sel.range <- "1-2"
      min_range = 1
      max_range = 2
      
    } else if (m_0_1_parameters[i, "MSE"] <= m_1_2_parameters[i, "MSE"]) {
      
      sel.parameters <- m_0_1_parameters[i, ]
      sel.values <- m_0_1_values[i, ]
      sel.derivatives <- m_0_1_derivatives[i, ]
      sel.range <- "0-1"
      min_range = 1e-5
      max_range = 1
      
    } else {
      
      sel.parameters <- m_1_2_parameters[i, ]
      sel.values <- m_1_2_values[i, ]
      sel.derivatives <- m_1_2_derivatives[i, ] 
      sel.range <- "1-2"
      min_range = 1
      max_range = 2
      
    }
    
    
    m_parameters <- rbind(m_parameters, setNames(sel.parameters, names(m_parameters)))
    m_values <- rbind(m_values, setNames(sel.values, names(m_values)))
    m_derivatives <- rbind(m_derivatives, setNames(sel.derivatives, names(m_derivatives)))
    range <- c(range, sel.range)
    
    x <- m_raw_values[i, 2:ncol(m_raw_values)]
    shifted_fit <- rescale_fun(min_range = min_range, 
                               max_range = max_range, 
                               x = x, 
                               y = as.numeric(sel.values[2:(tp+1)]))
    
    m_rescaled_values <- rbind(m_rescaled_values, setNames(c(m_raw_values[i, 1], shifted_fit), names(m_rescaled_values)))
    rescaled_b <- c(rescaled_b, rescale_fun(min_range = min_range, 
                                            max_range = max_range, 
                                            x = x, 
                                            y = as.numeric(t(sel.parameters[, 3]))))
  }
}  

m_parameters$range <- range
m_parameters$rescaled_b <- rescaled_b

colnames(m_rescaled_values) <- colnames(m_values)

# 3. save output
write.table(m_parameters, file = "parameters.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_values, file = "values.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_derivatives, file = "derivatives.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_rescaled_values, file = "rescaled.values.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
