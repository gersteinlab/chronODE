if (params.help) {
  log.info 'Usage: nextflow run chronode.nf [options]'
  log.info 'Parameters: '
  log.info ' --infile   Input matrix (should already be fully transformed)'
  log.info ' --size     Chunk size (affects speed but not final output)'
  log.info ' --out      Prefix for output files'
  log.info ' --dir      Output directory'
  log.info ' --timesfile  File of original timepoints' 
}

process fitting {
  input: 
  file x

  output:
  path('values.tsv')
  path('derivatives.tsv')
  path('parameters.tsv')
  path('rescaled.values.tsv')

  script:
  """
  source /vast/palmer/apps/avx2/software/miniconda/24.3.0/etc/profile.d/conda.sh
  conda activate "/gpfs/gibbs/pi/gerstein/bb926/conda_envs/jupiter"

  # get normalized matrices
  min_max.norm.R -i $x  
  
  # ODE fitting (range 0-1)
  ODE_fitting.py -t ${params.timesfile} -i minmax.0_1.tsv -v values.0_1.tsv -d derivatives.0_1.tsv -p parameters.0_1.tsv

  # ODE fitting (range 1-2)
  ODE_fitting.py -t ${params.timesfile} -i minmax.1_2.tsv -v values.1_2.tsv -d derivatives.1_2.tsv -p parameters.1_2.tsv

  # select best fitting
  select.fitting.R --p_0_1 parameters.0_1.tsv \
                   --v_0_1 values.0_1.tsv \
                   --d_0_1 derivatives.0_1.tsv \
                   --p_1_2 parameters.1_2.tsv \
                   --v_1_2 values.1_2.tsv \
                   --d_1_2 derivatives.1_2.tsv \
                   --raw_values $x
  """
}


workflow { 

    // split input matrix into chunks of a specified size
    split_ch = Channel.fromPath(params.infile) \
    | splitText(by: params.size, keepHeader: true, file: true)
    
    // run ODE fitting on each chunk
    (values_ch, derivatives_ch, parameters_ch, rescaled_values_ch) = fitting(split_ch)
    
    // merge values, parameters, derivatives, rescaled values across chunks
    collected_values = values_ch.collectFile(name: "${params.out}.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_derivatives = derivatives_ch.collectFile(name: "${params.out}.derivatives.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_parameters = parameters_ch.collectFile(name: "${params.out}.parameters.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_rescaled_values = rescaled_values_ch.collectFile(name: "${params.out}.rescaled.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
  
}

