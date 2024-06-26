# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based  pipeline for interpolating values and derivatives from time-series data, and machine learning models that predicted gene expression from linked chromatin features.

![](https://github.com/gersteinlab/chronODE/blob/main/figure1.png)

***

## Data interpolation and ODE fitting

#### Dependencies
```
argparse
numpy
pandas
scipy
```
#### Input requirements

Usage: `ODE_fitting_posvals.py [-h] [-i INPUTFILE] [-t TIMEPOINTS] [-T TIMECOURSE] [-g GROUP] [-r REGION] [-a ASSAY] [-v VALUESFILE] [-d DERIVSFILE]
                              [-p PARAMSFILE] [-f FOLDER]`
                                
Required arguments:
```
  -i INPUTFILE, --inputfile INPUTFILE
                        Input file with values at a handful of timepoints
  -v VALUESFILE, --valuesfile VALUESFILE
                        Values output file
  -d DERIVSFILE, --derivsfile DERIVSFILE
                        Derivatives output file
  -p PARAMSFILE, --paramsfile PARAMSFILE
                        Parameters output file
  -f FOLDER, --folder FOLDER
                        Base filepath for file names
```
Optional arguments:
```
  -h, --help            show this help message and exit

  -t TIMEPOINTS, --timepoints TIMEPOINTS
                        Number of timepoints to interpolate; defaults to 105
  -T TIMECOURSE, --timecourse TIMECOURSE
                        Set of initial timepoints to use, either mouse (default) or human
  -g GROUP, --group GROUP
                        expression pattern; defauts to 'unknown'
  -r REGION, --region REGION
                        Mouse brain region; defaults to 'unknown'
  -a ASSAY, --assay ASSAY
                        Assay suffix to put on column names; defaults to 'x'
```
#### File formats
The tab-separated input file needs one row per gene/regulatory element, one index column, and a column for each time point in the original data:  
```
cCRE_id	        E10.5	                E11.5	                E12.5	                E13.5	                E14.5	                E15.5	                E16.5	                PN
EM10D0043278	0.112209163706003	0.263615771996068	0.101058746899128	0.00109066358305005	-0.0733457337036519	-0.155362226485118	-0.130227596077729	-0.103008112890712
EM10D0046746	0.158520542990051	0.207235934578883	0.0882546996858981	0.0342933724159356	-0.0391589358988853	-0.136999257336671	-0.328550015608365	-0.181150203360637
```
The parameters output will be tab-separated and have a row for each element and have columns for the k and b parameters, mean squared error (MSE), fuction sign, the vertical move applied to the data, user-specified group, and user-specified region:  
```
cCRE_id	        k	                b	                MSE	                sign_func	TYPE	        MOVE	              group	  region
EM10D0043278	-0.053382675494526265	124805.313024191	0.02775699273405798	1.0	        upward	        0.9962693785260658    downreg	  forebrain
EM10D0046746	-0.8540763957605665	0.9343163272858892	0.011764908795083247	1.0	        original	0.0	              downreg	  forebrain
```
The derivatives and values output files will be tab-separated and have a row for each element and a column for each interpolated timepoint:
```
cCRE_id		10.5_oc			10.601_oc		10.702_oc		...	20.798_oc		20.899_oc		21.0_oc
EM10D0043278	-0.08727414570298016	-0.08680505073426538	-0.0863384770664966	...	-0.05036650875841838	-0.05009578784751738	-0.04982652204472581
EM10D0046746	-0.02097416583664084	-0.022752174901581188	-0.024670440353401078	...	-0.004305420230519206	-0.003953267439604704	-0.0036296506721193054
```
#### Examples
Run pipeline on RNAseq data from mouse forebrain:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/rna_example_input.tsv" \
      --timepoints 105 \
      --timecourse "mouse" \
      --region "forebrain" \
      --assay "RNAseq" \
      --valuesfile "example_dir/rna_vals.tsv" \
      --derivsfile "example_dir/rna_derivs.tsv" \
      --paramsfile "example_dir/rna_params.tsv"
Run pipeline on DNase-seq data from human brain decreasing cCREs:  

    python ODE_fitting_posvals.py \
      --inputfile "example_dir/human_dnase.tsv" \
      --timepoints 105 \
      --timecourse "human" \
      --group "decreasing" \
      --assay "DNase" \
      --valuesfile "example_dir/dnase_vals.tsv" \
      --derivsfile "example_dir/dnase_derivs.tsv" \
      --paramsfile "example_dir/dnase_params.tsv"

## Estimating switching time and classifying curves

#### Dependencies
```
optparse
```
#### Input requirements
`Usage: switching.time.labels.R [options]`  
Required arguments:
```
        -p PARAMSFILE, --paramsfile=PARAMSFILE
                Parameters file to use as input

        -v VALUESFILE, --valuesfile=VALUESFILE
                Values file to use as input

        -o OUTFILE, --outfile=OUTFILE
                Name of output file.
```
Optional arguments:
```

        -s START, --start=START
                Start time (in days), defaults to 10.5

        -e END, --end=END
                End time (in days), defaults to 21

        -g GROUP, --group=GROUP
                Pattern filter (increasing or decreasing), defaults to none

        -h, --help
                Show help message and exit
```
#### File formats
The input files should be the output values and parameters from a single chronODE run:
```
cCRE_id		10.5_oc			10.601_oc		10.702_oc		...	20.798_oc		20.899_oc		21.0_oc
EM10D0043278	-0.08727414570298016	-0.08680505073426538	-0.0863384770664966	...	-0.05036650875841838	-0.05009578784751738	-0.04982652204472581
EM10D0046746	-0.02097416583664084	-0.022752174901581188	-0.024670440353401078	...	-0.004305420230519206	-0.003953267439604704	-0.0036296506721193054
```
```
cCRE_id	        k	                b	                MSE	            	sign_func	TYPE	    	MOVE	              	group	  	region
EM10D0043278	-0.053382675494526265	124805.313024191	0.02775699273405798	1.0		upward	  	0.9962693785260658	downreg		forebrain
EM10D0046746	-0.8540763957605665	0.9343163272858892	0.011764908795083247	1.0		original	0.0	                downreg		forebrain
```
The output file is identical to the parameters input, but with columns added:
```
cCRE_id		k			b			MSE			sign_func	TYPE		MOVE			group		region		switching_time		saturation_time		minimum_time		label
EM10D0043278	-0.0533826754945263	124805.313024191	0.027756992734058	1		upward		0.996269378526066	decreasing	forebrain	-200.109828054027	-286.188688050577	490.027208851933	decelerator
EM10D0046746	-0.854076395760566	0.934316327285889	0.0117649087950832	1		original	0			decreasing	forebrain	14.6963669566587	9.31614585121227	57.8322757240436	switcher
```
#### Example
```
Rscript switching.time.labels.R \
    -p example_data/oc_params.tsv \
    -v example_data/oc_vals.tsv \
    -o example_data/oc_labeled.tsv \
    -s 10.5 \
    -e 21
```
## Predicting time-series gene expression velocities from chromatin velocity  of associated cCREs

### Random Forest

#### Dependencies
```
numpy
pandas
sklearn
```
#### Input requirements
`usage: randomforest.py [-h] [-i INPUTFILE] [-p PREDICTFILE] [-t TIMEPOINTS] [-m MODALITIES] [-s RANDOMSEED]`  
Required arguments:
```
  -i INPUTFILE, --inputfile INPUTFILE
                        Input file with paired RNA and chromatin feature values
  -p PREDICTFILE, --predictfile PREDICTFILE
                        File to write predictions to
```
Optional arguments:
```
  -h, --help            show help message and exit

  -t TIMEPOINTS, --timepoints TIMEPOINTS
                        Number of time points in the data; defaults to 105
  -m MODALITIES, --modalities MODALITIES
                        Number of input data modalities; defaults to 1
  -s RANDOMSEED, --randomseed RANDOMSEED
                        Random seed
```
#### File formats
The input file should be a tab-separated merge of the derivatives outputted by chronODE. Each row correspond to a linked gene-cCRE pair. The first two columns are expected to be indexes (i.e. the IDs of the gene-cCRE pairs). The next set of columns are the derivatives of chromatin accessibility (or other features) at each time point. The final set of columns are the derivatives of gene expression at each time point. Additions columns after the RNA columns will be ignored:
```
cCRE_id		gene_id			10.5_oc			10.601_oc		10.702_oc		...	20.798_oc		20.899_oc		21.0_oc			10.5_rna		10.601_rna		10.702_rna		...	20.899_rna		21.0_rna		correlation		k_rna			k_oc			region		group
EM10D1015674	ENSMUSG00000031138	0.0292021391502785	0.0310246038870716	0.0329472699456372	...	0.0148927361009084	0.0139712129144637	0.013104286814022	0.110106777956121	0.110665757107439	0.111212696967868	...	0.0887423120228607	0.0879574610481399	0.81216813024026	0.16939688531092	0.661612298960005	forebrain	late_upreg
EM10D1049003	ENSMUSG00000026085	0.284229691840175	0.281767561508326	0.278749667166499	...	0.0011924882022946	0.0011175292340872	0.0010472735390207	1.40716563847581	1.03119604214313	0.694939955232681	...	8.087689415355689e-16	0.0			0.489943751344613	5.37904223774371	0.644327248946338	forebrain	late_upreg
```
The output file will look like chronODE output, but only gene-cCRE pairs from the test set will be written to the predictions file. Note that column headers will also be replaced by generic integers:
```
combined_index					0			1			2			...	102			103			104
ENSMUSG00000039703-forebrain-EM10D3074493	-0.2625324998352258	-0.25478173365352913	-0.24735001435865087	...	-0.016483826507750404	-0.0160613323302765	-0.01564971378993606
ENSMUSG00000041895-midbrain-EM10D2386250	0.21748892698510805	0.16237279764122636	0.11835493752652373	...	0.0007211071451451441	0.0008984222916715194	0.0010708635950900517
```
#### Example
```
python randomforest.py \
  -i myfolder/rf_input.tsv \
  -p myfolder/rf_predictions.tsv \
  -s 42
 ```

### Neural Network

#### Dependencies
```
argparse
numpy
math
matplotlib
pandas
scipy
sklearn
torch
```
#### Input requirements
`usage: nn_model.py [-h] [-i INPUTFILE] [-p PREDICTFILE] [-S SUBSET] [-t TIMEPOINTS] [-n NEURALNETFILE] [-m MODALITIES] [-s SPLIT]
                   [-r RANDOMSEED]`
Required arguments:
```
  -i INPUTFILE, --inputfile INPUTFILE
                        Input file with paired RNA and OC values
  -p PREDICTFILE, --predictfile PREDICTFILE
                        File to write predictions to
  -S SUBSET, --subset SUBSET
                        Subset to model (choose "pos" or "neg")
```
Optional arguments:
```
  -h, --help            show this help message and exit
  -t TIMEPOINTS, --timepoints TIMEPOINTS
                        Number of time points in the data, defaults to 105
  -n NEURALNETFILE, --neuralnetfile NEURALNETFILE
                        Filepath to save model object, if used must end in .pth
  -m MODALITIES, --modalities MODALITIES
                        Number of input data modalities, defaults to 1
  -s SPLIT, --split SPLIT
                        Fraction of input to use as test set, defaults to 0.2
  -r RANDOMSEED, --randomseed RANDOMSEED
                        Random seed, defaults to 1941
```
#### File formats 
The input should be a table with columns for the IDs of a linked cCRE-gene pair, and columns for the cCRE's chromatin accessibility at each time point, followed by the gene's RNA epression levels at each time point:
```
cCRE_id 	gene_id 		10.5_oc 		10.601_oc       	10.702_oc       	...     20.798_oc       	20.899_oc       	21.0_oc 		10.5_rna        	10.601_rna      	10.702_rna      	...      20.798_rna      	20.899_rna		21.0_rna 		correlation     	k_rna   		k_oc    		region  	group
EM10D1015674    ENSMUSG00000031138      0.0292021391502785      0.0310246038870716      0.0329472699456372      ...     0.0148927361009084      0.0139712129144637      0.013104286814022       0.110106777956121       0.110665757107439       0.111212696967868       ...      0.0895245630016425     0.0887423120228607      0.0879574610481399      0.81216813024026        0.16939688531092  	0.661612298960005       forebrain       late_upreg
EM10D1048963    ENSMUSG00000026087      0.211646182284525       0.210998913384704       0.210139371063893       ...	0.0073492231516574      0.0070292450824912      0.0067229702510647      0.0182494031757812      0.0183168373878034      0.0183845207730992      ...      0.0265848293785538     0.0266830631330007      0.0267816598588292      -0.954902750003024      0.0365329628689712 	0.448571983858911       forebrain       late_upreg
EM10D1049003    ENSMUSG00000026085      0.284229691840175       0.281767561508326       0.278749667166499       ...     0.0011924882022946      0.0011175292340872      0.0010472735390207      1.40716563847581        1.03119604214313        0.694939955232681       ...      0       		8.08768941535569e-16    0       		0.489943751344613  	5.37904223774371        0.644327248946338       forebrain       late_upreg
```
The main output of the program is a tab-separated table of predicted RNA expression derivatives with a row for each gene in the test set and a column for each time point:
```
0		1		2		...	102		103		104
-0.088076845	-0.0883278	-0.08855562	...	0.0065826476	0.00670442	0.006819457
0.0044836476	0.003699176	0.0027732253	...	0.008761667	0.00876233	0.008762874
-0.0012679771	-0.0022765324	-0.003380455	...	0.008090496	0.008154936	0.008213259
```
If the -n option is used, the model state will be saved as a `.pth` object that is loadable using `torch`.
#### Example
```
python nn_model.py \
    --inputfile "myfolder/nn_poscorr_input.tsv" \
    -p "myfolder/nn_pos_results.tsv" \
    -S "pos"
```
## Euler Method Reconstruction

#### Dependencies
```
argparse
numpy
pandas
```
#### Input requirements
`usage: euler_method.py [-h] [-v VALUESFILE] [-p PREDICTFILE] [-o OUTFILE] [-t TIMEPOINTS] [-m MODE] [-g GENEID]`

Required arguments:
```
  -v VALUESFILE, --valuesfile VALUESFILE
                        Input file with actual values
  -p PREDICTFILE, --predictfile PREDICTFILE
                        Input file with predicted derivatives
  -o OUTFILE, --outfile OUTFILE
                        Output filepath
```
Optional arguments:
```
  -h, --help            show this help message and exit
  -t TIMEPOINTS, --timepoints TIMEPOINTS
                        Number of time points in the data, defaults to 105
  -m MODE, --mode MODE  Running mode, defaults to "all"; use "gene" for single gene mode
  -g GENEID, --geneid GENEID
                        Gene ID to use with single gene mode
```
#### File formats
The file of "true values" as fitted using an ODE needs a column for gene ID and a column for each timepoint:
```
gene_id 		10.5_rna        	10.601_rna      	10.702_rna      	...      20.798_rna      	20.899_rna      	21.0_rna
ENSMUSG00000000001      2.0000001334444906      1.9948365616249952      1.9896002014520098      ...      1.1028228229287411     1.0917646243698502      1.0807178757568114
ENSMUSG00000000028      2.000000261578661       1.9921621589473684      1.9842687414094864      ...      1.0321240089552974     1.0226041127895504      1.013118957079147
ENSMUSG00000000031      1.0000009999989565      0.9771926384317687      0.954904497183259       ...      0.09504614914040505    0.09287830237220789     0.0907599006284299
```
The file of modeled derivatives predicted by either the random forest or the neural net needs a column called combined_index, and a column for each timepoint:
```
combined_index 					0       		1       		2       		...     102     		103     		104
ENSMUSG00000090264-forebrain-EM10D2584003       -0.021443956144045807   -0.020675566736395985   -0.019959201046996285   ...     0.0423777198170062      0.0425623879262709      0.042741987452785006
ENSMUSG00000072620-midbrain-EM10D1246505        0.3550414284830037      0.2634913635363989      0.21427529099667958    	...   	0.00681893514174904     0.006657113064798876    0.006500109235992046
ENSMUSG00000069743-hindbrain-EM10D2550322       -0.15298214703368593    -0.15004397764262156    -0.14715171745341837    ...   	-0.005825759676286866   -0.005518792516383729   -0.0052189183160649515
```
The output format depends on the running mode selected. Single-gene mode produces a text file with the reconstructed gene expression value at each timepoint on its own line:
```
2.0000009999597714
1.9864202534124273
1.970736494777169
1.9528113120211568
```
All mode produces a tab-separated table with a column for each timepoint and a row for each gene shared that was found in both input sets:
```
gene_id			0			1			2			...	102			103			104
ENSMUSG00000001039	2.0000009999597714	1.9864202534124273	1.970736494777169	...	1.8114385194575986	1.8130753497051955	1.8147339436463965
ENSMUSG00000001366	1.0000005112607953	1.0131987265407805	1.0262204390773904	...	1.2356192369446632	1.2340631524997345	1.2326833794704843
ENSMUSG00000001472	0.9998183996685308	1.0304834388749835	1.0513223972192403	...	1.8276526132574102	1.827656937757639	1.8276608899911557
```
#### Example
Using single-gene mode:
```
python euler_method.py \
	--valuesfile "example_dir/rna_vals.tsv" \
	--predictfile "myfolder/rf_predictions.tsv" \
	--outfile "myfolder/rf_euler.txt" \
	-m "gene"
	--gene TK
```
Using all mode:
```
python euler_method.py \
	--valuesfile "example_dir/rna_vals.tsv" \
	--predictfile "myfolder/rf_predictions.tsv" \
	--outfile "myfolder/rf_euler.tsv" \
	--mode "all"
```
## License
Copyright 2024 Yale University.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


