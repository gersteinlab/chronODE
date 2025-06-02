# The chronODE framework for modelling multi-omic time series with ordinary differential equations and machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based framework for modelling logistic and peak-like time-series multi-omic signals, and machine learning models that predicted gene expression over time from linked chromatin features.

![](https://github.com/gersteinlab/chronODE/blob/main/figure1.png)

***

## Monotonic fitting

#### Dependencies

The chronODE pipeline for monotonic fitting is written in [Nextflow](https://www.nextflow.io/) DSL2.
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
The tab-separated main input file needs one row per gene/regulatory element, one index column, and a column for each time point in the original data. Below, we provide an example for two candidate cis-regulatory elements (cCREs):
```
cCRE_id	        E10.5	                E11.5	                E12.5	                E13.5	                E14.5	                E15.5	                E16.5	                PN
EM10D0144246	1.92265260038112	1.89139781357557	1.90998930285207	2.24724966250699	2.44195099910917	2.56235923617873	2.40102884091981	2.5533492607186
EM10D1047237	1.44726021325062	1.39525207794537	1.41927175130354	1.3478911616366		1.16733471629674	1.09028685205392	1.08475277892454	1.07986151310049
```

The time course must be specified using a .csv file listing all time points in numeric form on one line. Please consider that in our example we have numerically encoded the first post-natal day (e.g., PN) as day 21.
```
10.5,11.5,12.5,13.5,14.5,15.5,16.5,21
```

The parameters output will be tab-separated and have a row for each gene/regulatory element:
```
cCRE_id       k                  b_starred          MSE                  a                 b                 R_min  R_max  z_min             z_max             z_start           kinetic_class  switching_time    saturation_time   minimum_time
EM10D0144246  1.02066198626769   0.972883868461625  0.0162867863978315   1.89139110389425  2.54416517604565  1e-05  1      1.89139781357557  2.56235923617873  1.92265260038112  switcher       13.429223902119   17.9313214740467  -22.6663317113752
EM10D1047237  -1.99779861566745  1.00103040140379   0.00245095219378508  1.07985783907675  1.44763878517272  1e-05  1      1.07986151310049  1.44726021325062  1.44726021325062  switcher       13.9426932656304  11.6426016477376  32.3836718502493
```
This file contains the following information:
* `k`, `a`, and `b`: the three logistic ODE parameters in the real range of the data.
* `b_starred`: the upper asymptote of the logistic curve in the normalized range of the data (simplified form of the logistic ODE).
* `R_min`, `R_max`, `z_min`, and `z_max`: scaling factors used in the data normalization step (see Methods section of the manuscript "Data normalization").
* `z_start`: the gene/element's signal at the first experimental time point in the real range of the data. 
* `MSE`: the Mean Squared Error of the fit in the normalized range of the data. 
* `kinetic_class`, `switching_time`, `saturation_time`, and `minimum_time`: kinetic characterization of the element (see Methods, section `2.4 Kinetic Classification`, and Supplementary Note 1).

The derivatives, fitted values, and restored values output files will also be tab-separated and have a row for each gene/element and a column for each timepoint. Please see below an example of the output file for the `restored.values.tsv` output file:
```
cCRE_id       E10.5             E11.5             E12.5             E13.5             E14.5             E15.5             E16.5             PN
EM10D0144246  1.92265260038112  2.01453994115168  2.22955836852258  2.43480368195092  2.51693257926945  2.53807509249146  2.54283806152758  2.54387763550249
EM10D1047237  1.44726021325062  1.44020555489426  1.34014945343893  1.11954446308323  1.08206670365814  1.07996880889247  1.07986338380165  1.07985811604938
```
#### Example
```
nextflow run chronode.nf \
  --infile example_data/input/input.matrix.tsv \
  --out example.output \
  --size 3 \
  --dir example_data/output \
  --timesfile example_data/input/mouse.timecourse.csv
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


