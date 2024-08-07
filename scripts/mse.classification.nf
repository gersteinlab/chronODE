if (params.help) {
  log.info 'Usage: nextflow run chronode.nf [options]'
  log.info 'Parameters: '
  log.info ' --orig            Matrix of original values'
  log.info ' --fitted          Matrix of fitted values (0-1 or 1-2 range)'
  log.info ' --rescfitted      Matrix of fitted rescaled values (original range)'
  log.info ' --deriv           Matrix of derivatives'
  log.info ' --par             Matrix of parameters'
  log.info ' --tokeep          List of reproducibly fitted elements'
  log.info ' --size            Chunk size (affects speed but not final output)'
  log.info ' --tstart          Initial time point of the time course'
  log.info ' --tend            Last time point of the time course'
  log.info ' --out             Prefix for output files'
  log.info ' --dir             Output directory'
}

process filtering {

  input: 
  tuple file(orig), file(fitted), file(rescfitted), file(deriv), file(par), file(tokeep)

  output:

  path("filtered.orig.tsv") 
  path("filtered.fitted.tsv") 
  path("filtered.rescfitted.tsv")
  path("filtered.deriv.tsv")
  
  tuple path("filtered.orig.tsv"), path("filtered.fitted.tsv"), path("filtered.rescfitted.tsv"), path("filtered.par.tsv") 

  script:
  """
  source /vast/palmer/apps/avx2/software/miniconda/24.3.0/etc/profile.d/conda.sh
  conda activate "/gpfs/gibbs/pi/gerstein/bb926/conda_envs/jupiter"

  filter.R -o $orig -f $fitted -r $rescfitted -d $deriv -p $par -k $tokeep
  """

}

process mse_classification {
  
  input:
  tuple file(f_orig), file(f_fitted), file(f_rescfitted), file(f_par)

  output:
  path("mseClass.param.tsv")

  script:
  """
  source /vast/palmer/apps/avx2/software/miniconda/24.3.0/etc/profile.d/conda.sh
  conda activate "/gpfs/gibbs/pi/gerstein/bb926/conda_envs/jupiter"

  # compute MSE permutation p-value
  mse.perm.pval.py -o $f_orig -r $f_rescfitted -p $f_par
  
  # add FDR and classification of elements (switchers / accelerators / decelerators)
  switching.time.labels.R -p mse.param.tsv -f $f_fitted -o mseClass.param.tsv -s ${params.tstart} -e ${params.tend}
  """

}

workflow { 

    // create input channels  
    input_ch = Channel.of([params.orig, params.fitted, params.rescfitted, params.deriv, params.par, params.tokeep]) \
       | map { orig, fitted, rescfitted, deriv, par, tokeep -> tuple(
        file(orig),
        file(fitted),
        file(rescfitted),
        file(deriv),
        file(par),
        file(tokeep)) } 


    // filter elements based on reproducibility list & split files into chunks
    (orig_ch, fitted_ch, rescfitted_ch, deriv_ch, mse_ch) = filtering(input_ch) 


    input_ch = mse_ch.flatten() \
     | splitText(by: params.size, keepHeader: true, file: true) \
     | flatten \
     | map { it ->
         def it2 = it.toString().tokenize('/').last().replace('.tsv', '')
         def f = it2.tokenize(".").last()
         [f, it] } \
     | groupTuple(by: [0]) \
     | map { groupTuple ->
         def secondElement = groupTuple[1]
         return secondElement } 


    // run MSE p-value calculation and classification on each chunk
    param_ch  = mse_classification(input_ch)


    // merge values, parameters, derivatives, rescaled values across chunks
    collected_orig = orig_ch.collectFile(name: "${params.out}.original.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_fitted = fitted_ch.collectFile(name: "${params.out}.fitted.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_rescaled_values = rescfitted_ch.collectFile(name: "${params.out}.rescaled.values.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_derivatives = deriv_ch.collectFile(name: "${params.out}.derivatives.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")
    collected_parameters = param_ch.collectFile(name: "${params.out}.parameters.tsv", keepHeader: true, skip: 1, storeDir: "${params.dir}")

}