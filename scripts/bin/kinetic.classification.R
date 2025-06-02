#!/usr/bin/env Rscript


#***************************
# author: Beatrice Borsari |
#***************************

# This script selects the best fitting between 
# the fittings performed in the 0-1 and 1-2 range
#
# For the selected fitting, it then computes 
# a and b in the original range of the data
# and it uses these parameters to restore the fitted curve
# to the original range of the data.
#
# Finally, it performs a kinetic classification of the elements
# and computes switching time, minimum time, and saturation time.



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
  make_option( c ( "--v_0_1" ), help = "The 0-1 fitted values matrix (y*_fitted)."),
  make_option( c ( "--v_1_2" ), help = "The 1-2 fitted values matrix (y*_fitted)."),
  make_option( c ( "--d_0_1" ), help = "The 0-1 derivatives matrix."),
  make_option( c ( "--d_1_2" ), help = "The 1-2 derivatives matrix."),
  make_option( c ( "--time_points" ), help = "The list of time points used."),
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

# this function is used to rescale the parameters
# a and b back to the original range of the data
# (this step is described in Methods section 2.1
# "Restoring the fitted curve to the original range of the data.")

rescale_fun <- function(f_R_min, f_R_max, f_z_min, f_z_max, par_fitted) {
  
  par_restored = (((par_fitted - f_R_min)*(f_z_max - f_z_min))/(f_R_max - f_R_min)) + f_z_min
  
  return(par_restored)
  
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
#             raw_values = "forebrain.dynamic.global.shift.6.cCREs.2.tsv",
#             time_points = "/gpfs/gibbs/pi/gerstein/bb.shared.projects/brain-comp-dev/analyses/mouse.NatComm.rev1/ODE/fb.mb.hb/ODE.fitting/mouse.timecourse.csv")

#------------------------
# 1. read in input data |
#------------------------
# 1.1. parameters in range 0-1 and 1-2
m_0_1_parameters <- read.table(opt$p_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_parameters <- read.table(opt$p_1_2, header = T, sep="\t", stringsAsFactors = F)

# # only for debugging
# colnames(m_0_1_parameters)[3] <- "b_starred"
# colnames(m_1_2_parameters)[3] <- "b_starred"


# 1.2. fitted values in range 0-1 and 1-2 (i.e., this is y*_fitted or y_starred_fitted)
m_0_1_values <- read.table(opt$v_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_values <- read.table(opt$v_1_2, header = T, sep="\t", stringsAsFactors = F)


# 1.3. derivatives (these are the derivatives of y*_fitted)
m_0_1_derivatives <- read.table(opt$d_0_1, header = T, sep="\t", stringsAsFactors = F)
m_1_2_derivatives <- read.table(opt$d_1_2, header = T, sep="\t", stringsAsFactors = F)


# 1.4. raw values (these are the z values)
m_raw_values <- read.table(opt$raw_values, header = T, sep="\t", stringsAsFactors = F)


# 1.5. time points
timepoints <- read.table(opt$time_points, h=F, sep = ",")
t_start <- timepoints[1, 1]
t_end <- timepoints[1, ncol(timepoints)]

# get number of time points from fitted values
tp_len = ncol(m_0_1_values) - 1


#------------------------------------------------------------------
# 2. select best fitting for each element                         |
# note: the parameters df at the beginning contains only 3 values |
# k, b* (b_starred), and MSE                                      |
# at the end of the script, we will add the following info:       |
# - a & b                                                         |
# - R_min & R_max                                                 |
# - z_min & z_max                                                 |
# - z_start                                                       |
# - switching_time, minimum_time, saturation_time                 |
# - kinetic classification                                        |
#------------------------------------------------------------------

# 2.1. define objects to store output values
m_parameters <- data.frame(matrix(ncol=ncol(m_0_1_parameters),nrow=0, dimnames=list(NULL, colnames(m_0_1_parameters))), stringsAsFactors = F)
m_fitted_values <- data.frame(matrix(ncol=ncol(m_0_1_values),nrow=0, dimnames=list(NULL, colnames(m_0_1_values))), stringsAsFactors = F)
m_derivatives <- data.frame(matrix(ncol=ncol(m_0_1_derivatives),nrow=0, dimnames=list(NULL, colnames(m_0_1_derivatives))), stringsAsFactors = F)
m_restored_values <- data.frame(matrix(ncol=ncol(m_0_1_values),nrow=0, dimnames=list(NULL, colnames(m_0_1_values))), stringsAsFactors = F)
a <- c()
b <- c()
R_min <- c()
R_max <- c()
z_min <- c()
z_max <- c()
z_start <- c()
minimum_time <- c()
saturation_time <- c()
switching_time <- c()
kinetic_class <- c()

# 2.2. perform calculations
for ( i in 1:nrow(m_1_2_parameters) ) {
  
  # 2.2.1. check if the element has an available fitting either in the 0-1 or 1-2 range
  if (is.na(m_1_2_parameters[i, "MSE"]) & (is.na(m_0_1_parameters[i, "MSE"]))) {
    
    # 2.2.1.1. if there's no available fitting, we will add NA for this element
    m_parameters <- rbind(m_parameters, data.frame(t(setNames(c(m_1_2_parameters[i, 1], rep(NA, 3)), names(m_parameters)))))
    m_fitted_values <- rbind(m_fitted_values, data.frame(t(setNames(c(m_1_2_values[i, 1], rep(NA, tp_len)), names(m_fitted_values)))))
    m_derivatives <- rbind(m_derivatives, data.frame(t(setNames(c(m_1_2_derivatives[i, 1], rep(NA, tp_len)), names(m_derivatives)))))
    m_restored_values <- rbind(m_restored_values, data.frame(t(setNames(c(m_1_2_values[i, 1], rep(NA, tp_len)), names(m_restored_values)))))
    
    a <- c(a, NA)
    b <- c(b, NA)
    R_min <- c(R_min, NA)
    R_max <- c(R_max, NA)
    z_min <- c(z_min, NA)
    z_max <- c(z_max, NA)
    z_start <- c(z_start, NA)
    minimum_time <- c(minimum_time, NA)
    saturation_time <- c(saturation_time, NA)
    switching_time <- c(switching_time, NA)
    kinetic_class <- c(kinetic_class, NA)

    # 2.2.1.2. otherwise, we pick the fitting that is available and that has the lowest MSE
  } else {
    
    ## fitting available only in the 0-1 range
    if (is.na(m_1_2_parameters[i, "MSE"]) & !(is.na(m_0_1_parameters[i, "MSE"]))) {
      
      sel.parameters <- m_0_1_parameters[i, ]
      sel.fitted.values <- m_0_1_values[i, ]
      sel.derivatives <- m_0_1_derivatives[i, ]
      
      i_R_min = 1e-5
      i_R_max = 1
      
    ## fitting available only in the 1-2 range
    } else if (is.na(m_0_1_parameters[i, "MSE"]) & !(is.na(m_1_2_parameters[i, "MSE"]))) {
      
      sel.parameters <- m_1_2_parameters[i, ]
      sel.fitted.values <- m_1_2_values[i, ]
      sel.derivatives <- m_1_2_derivatives[i, ] 

      i_R_min = 1
      i_R_max = 2
      
    ## fitting in the 0-1 range has lower MSE than fitting in the 1-2 range
    } else if (m_0_1_parameters[i, "MSE"] <= m_1_2_parameters[i, "MSE"]) {
      
      sel.parameters <- m_0_1_parameters[i, ]
      sel.fitted.values <- m_0_1_values[i, ]
      sel.derivatives <- m_0_1_derivatives[i, ]

      i_R_min = 1e-5
      i_R_max = 1
      
    ## fitting in the 1-2 range has lower MSE than fitting in the 0-1 range
    } else {
      
      sel.parameters <- m_1_2_parameters[i, ]
      sel.fitted.values <- m_1_2_values[i, ]
      sel.derivatives <- m_1_2_derivatives[i, ] 

      i_R_min = 1
      i_R_max = 2
      
    }
    
    
    # 2.2.1.3. record the parameters, fitted values, and derivatives in the selected range
    m_parameters <- rbind(m_parameters, setNames(sel.parameters, names(m_parameters)))
    m_fitted_values <- rbind(m_fitted_values, setNames(sel.fitted.values, names(m_fitted_values)))
    m_derivatives <- rbind(m_derivatives, setNames(sel.derivatives, names(m_derivatives)))

    
    # 2.2.1.4. record selected range [R_min-Rmax]
    R_min <- c(R_min, i_R_min)
    R_max <- c(R_max, i_R_max)
    
    
    # 2.2.1.5. record z_start, z_min and z_max
    z <- as.numeric(m_raw_values[i, 2:ncol(m_raw_values)])
    i_z_min <- min(z)
    i_z_max <- max(z)
    i_z_start <- z[1]
    
    z_min <- c(z_min, i_z_min)
    z_max <- c(z_max, i_z_max)
    z_start <- c(z_start, i_z_start)

    
    # 2.2.1.6. compute b in the real range of the data
    i_b <- rescale_fun(f_R_min = i_R_min, 
                       f_R_max = i_R_max, 
                       f_z_min = i_z_min,
                       f_z_max = i_z_max,
                       par_fitted = as.numeric(t(sel.parameters[, "b_starred"])))
    b <- c(b, i_b)
    
    
    # 2.2.1.7. compute a in the real range of the data
    i_a <- rescale_fun(f_R_min = i_R_min, 
                       f_R_max = i_R_max, 
                       f_z_min = i_z_min,
                       f_z_max = i_z_max,
                       par_fitted = 0)
    a <- c(a, i_a)
    
    
    # 2.2.1.8. restore the fitted curve back to the real range of the data
    i_L = i_b - i_a
    i_k = as.numeric(t(sel.parameters[, "k"]))
    i_c = log(abs((i_z_start - i_a) / (i_L - (i_z_start - i_a)))) - i_k * t_start
    i_C = exp(i_c)
    tp_out = seq(t_start, t_end, length.out = tp_len)
    
    i_ID = sel.parameters[1, 1]
    
    i_restored_values = c(i_ID, 
                         (((i_L*i_C*exp(i_k * tp_out)) / (1 + i_C*exp(i_k * tp_out))) + i_a))
    
    m_restored_values <- rbind(m_restored_values, setNames(i_restored_values, names(m_restored_values)))
    
    
    # 2.2.1.9. get b_starred and y_starred_start
    i_b_starred <- as.numeric(t(sel.parameters[, "b_starred"]))
    i_y_starred_start <- as.numeric(t(sel.fitted.values[, 2]))
    
    
    # 2.2.1.10. compute C_starred
    i_c_starred = log(i_y_starred_start) - log(1 - (i_y_starred_start/i_b_starred)) - (i_k * t_start)
    i_C_starred = exp(i_c_starred)
    
    
    # 2.2.1.11. compute t_switch, t_min, t_sat
    i_switching_time <- (log(i_b_starred/i_C_starred)/i_k)
    switching_time <- c(switching_time, i_switching_time)
    
    minimum_time <- c(minimum_time, 
                      (log(((10^-16)*(i_b_starred^2))/((i_b_starred*i_C_starred)-((10^-16)*i_b_starred*i_C_starred)))/i_k))
    
    saturation_time <- c(saturation_time, 
                         (log((0.99*(i_b_starred^2))/((i_b_starred*i_C_starred)-(0.99*i_b_starred*i_C_starred)))/i_k))
    
    
    # 2.2.1.12. get kinetic classification
    
    if (i_switching_time > t_end) 
      
      {kinetic_class <- c(kinetic_class, "accelerator")} 
    
    else if (i_switching_time < t_start) 
      
      {kinetic_class <- c(kinetic_class, "decelerator")} 
    
    else 
      
      {kinetic_class <- c(kinetic_class, "switcher")}
    
  }
}  


# 3. add newly-computed parameters to the parameters df
m_parameters$a <- a
m_parameters$b <- b
m_parameters$R_min <- R_min
m_parameters$R_max <- R_max
m_parameters$z_min <- z_min
m_parameters$z_max <- z_max
m_parameters$z_start <- z_start
m_parameters$kinetic_class <- kinetic_class
m_parameters$switching_time <- switching_time
m_parameters$saturation_time <- saturation_time
m_parameters$minimum_time <- minimum_time


# 4. update colnames of m_restored_values
colnames(m_restored_values) <- colnames(m_fitted_values)


# 5. save output
write.table(m_parameters, file = "parameters.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_fitted_values, file = "fitted.values.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_derivatives, file = "derivatives.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
write.table(m_restored_values, file = "restored.values.tsv",
            row.names = F, col.names = T, sep="\t", quote = F)
