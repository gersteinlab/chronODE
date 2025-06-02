# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based  pipeline for interpolating values and derivatives from time-series data, and machine learning models that predicted gene expression from linked chromatin features.

![](https://github.com/gersteinlab/chronODE/blob/main/figure1.png)

***

## Monotonic fitting

#### Dependencies

The chronode pipeline for monotonic fitting is written in [Nextflow](https://www.nextflow.io/) DSL2.
The list of dependencies is provided in the `chronode.yml` file (see above). 
To create the corresponding conda environment, please run: 

```
conda env create -f chronode.yml
```

#### Input requirements

Usage: `nextflow run chronode.nf [options]`

Parameters: 
```
 --infile     Input matrix 
 --size       Chunk size (affects speed but not final output). We recommend ~10% of the file length, depending on available resources.
 --out        Prefix for output files
 --dir        Output directory
 --timesfile  File of time-points
```

Nextflow offers a number of other optional parameters, including `--help` and `-with-trace` that may be useful for debugging if errors occur.
#### File formats
The tab-separated main input file needs one row per gene/regulatory element, one index column, and a column for each time point in the original data:  
```
cCRE_id	        E10.5	                E11.5	                E12.5	                E13.5	                E14.5	                E15.5	                E16.5	                PN
EM10D0043278	0.112209163706003	0.263615771996068	0.101058746899128	0.00109066358305005	-0.0733457337036519	-0.155362226485118	-0.130227596077729	-0.103008112890712
EM10D0046746	0.158520542990051	0.207235934578883	0.0882546996858981	0.0342933724159356	-0.0391589358988853	-0.136999257336671	-0.328550015608365	-0.181150203360637
```
The time course must be specified using a .csv file listing all time points in numeric form on one line:
```
10.5,11.5,12.5,13.5,14.5,15.5,16.5,21
```

The parameters output will be tab-separated and have a row for each element and have columns for the k and b parameters fitted to the data, mean squared error (MSE), the range that the data were shifted into for the fitting, and a rescaled value of b based on the original range.
```
cCRE_id		k			b			MSE			range	rescaled_b
EM10D2246738	0.120000858009469	239089.428129566	0.0557348329407984	0-1	127612.186705749
EM10D2246742	1.03068354441858	0.886385177421958	0.0237505021001379	0-1	2.88606678854373
```
The derivatives, fitted values, and rescaled fitted values output files will be tab-separated and have a row for each element and a column for each interpolated timepoint:
```
cCRE_id		10.5			11.5			12.5			13.5			14.5			15.5			16.5			21.0
EM10D1138540	1.51000458715981	1.43871659514363	1.34337943231759	1.23596323549907	1.13599611554697	1.05821095446409	1.00569653913472	0.937419990935623
EM10D1138541	1.70323032149406	1.40624378601762	1.20345995048573	1.06499817485729	0.970455811708974	0.905901834291343	0.861824070290061	0.783982973335906
```
#### Example
```
nextflow run chronode.nf \
  --infile example_input.tsv \
  --out forebrain_test \
  --size 3 \
  --dir example_output/ \
  --timesfile mouse.timecourse.csv
```

## Piecewise fitting
```
python piecewise.fitting.py
```

## Temporal prediction of gene expression

### bidirectional Recurrent Neural Network (biRNN)
```
python chronODE_biRNN_model.py
```
Inside ```models/biRNN```, there is a separate file for each of the four trained models (Enhancer Monopattern, Enhancer Polypattern, Silencer Monopattern, Silencer Polypattern).


## License
Copyright 2024 Yale University.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


