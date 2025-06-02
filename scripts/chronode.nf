/*
 * Copyright (c) 2025, Yale University
 *
 * This file is part of 'chronode-nf':
 * A Nextflow pipeline for modelling of multi-omic time-series.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */



// Print Usage

if (params.help) {
  log.info 'Usage: nextflow run chronode.nf [options]'
  log.info 'Parameters: '
  log.info ' --infile   Input matrix'
  log.info ' --size     Chunk size (affects speed but not final output)'
  log.info ' --out      Prefix for output files'
  log.info ' --dir      Output directory'
  log.info ' --timesfile  File of original timepoints' 
}



// Define main process

process fitting {
  input: 
  file x

  output:
  path('fitted.values.tsv')
  path('derivatives.tsv')
  path('parameters.tsv')
  path('restored.values.tsv')

  script:
  """

  # get normalized matrices
  min_max.norm.R -i $x  
  
  # ODE fitting (range 0-1)
  ODE_fitting.py -t ${params.timesfile} -i minmax.0_1.tsv -v values.0_1.tsv -d derivatives.0_1.tsv -p parameters.0_1.tsv

  # ODE fitting (range 1-2)
  ODE_fitting.py -t ${params.timesfile} -i minmax.1_2.tsv -v values.1_2.tsv -d derivatives.1_2.tsv -p parameters.1_2.tsv

  # select best fitting, compute a and b parameters in the original range, restored fitted curve to the original range,
  # and perform kinetic classification 
  kinetic.classification.R --p_0_1 parameters.0_1.tsv \
			   --v_0_1 values.0_1.tsv \
			   --d_0_1 derivatives.0_1.tsv \
			   --p_1_2 parameters.1_2.tsv \
			   --v_1_2 values.1_2.tsv \
			   --d_1_2 derivatives.1_2.tsv \
			   --time_points ${params.timesfile} \
			   --raw_values $x
  """
}


workflow { 

    // split input matrix into chunks of a specified size
    split_ch = Channel.fromPath(params.infile) \
    | splitText(by: params.size, keepHeader: true, file: true)
    
    // run ODE fitting on each chunk
    (fitted_values_ch, derivatives_ch, parameters_ch, restored_values_ch) = fitting(split_ch)
    
    // merge fitted values, parameters, derivatives, and restored values across chunks
    collected_fitted_values = fitted_values_ch.collectFile(name: "${params.out}.fitted.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_derivatives = derivatives_ch.collectFile(name: "${params.out}.derivatives.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_parameters = parameters_ch.collectFile(name: "${params.out}.parameters.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_restored_values = restored_values_ch.collectFile(name: "${params.out}.restored.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
  
}

