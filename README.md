# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based  pipeline for interpolating values and derivatives from time-series data, and machine learning models that predicted gene expression from linked chromatin features.

![](https://github.com/gersteinlab/chronODE/blob/main/figure1.png)

***

## Data interpolation and ODE fitting

### Dependencies
```
argparse
numpy
pandas
scipy
```
### Input requirements

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
### File formats
The tab-separated input file needs one row per gene/regulatory element, one index column, and a column for each time point in the original data:  
```
cCRE_id	        E10.5	                E11.5	                E12.5	                E13.5	                E14.5	                E15.5	                E16.5	                PN
EM10D0043278	0.112209163706003	0.263615771996068	0.101058746899128	0.00109066358305005	-0.0733457337036519	-0.155362226485118	-0.130227596077729	-0.103008112890712
EM10D0046746	0.158520542990051	0.207235934578883	0.0882546996858981	0.0342933724159356	-0.0391589358988853	-0.136999257336671	-0.328550015608365	-0.181150203360637
```
The parameters output will be tab-separated and have a row for each element and have columns for the k and b parameters, mean squared error (MSE), fuction sign, the vertical move applied to the data, user-specified group, and user-specified region:  
```
cCRE_id	        k	                    b	                  MSE	            sign_func	TYPE	    MOVE	              group	  region
EM10D0043278	-0.053382675494526265	124805.313024191	 0.02775699273405798	  1.0	upward	  0.9962693785260658	downreg	forebrain
EM10D0046746	-0.8540763957605665	  0.9343163272858892	0.011764908795083247	1.0	original	0.0	                downreg	forebrain
```
The derivatives and values output files will be tab-separated and have a row for each element and a column for each interpolated timepoint:
```
cCRE_id	10.5_oc	10.601_oc	10.702_oc	10.803_oc	10.904_oc	11.005_oc	11.106_oc	11.207_oc	11.308_oc	11.409_oc	11.51_oc	11.611_oc	11.712_oc	11.812_oc	11.913_oc	12.014_oc	12.115_oc	12.216_oc	12.317_oc	12.418_oc	12.519_oc	12.62_oc	12.721_oc	12.822_oc	12.923_oc	13.024_oc	13.125_oc	13.226_oc	13.327_oc	13.428_oc	13.529_oc	13.63_oc	13.731_oc	13.832_oc	13.933_oc	14.034_oc	14.135_oc	14.236_oc	14.337_oc	14.438_oc	14.538_oc	14.639_oc	14.74_oc	14.841_oc	14.942_oc	15.043_oc	15.144_oc	15.245_oc	15.346_oc	15.447_oc	15.548_oc	15.649_oc	15.75_oc	15.851_oc	15.952_oc	16.053_oc	16.154_oc	16.255_oc	16.356_oc	16.457_oc	16.558_oc	16.659_oc	16.76_oc	16.861_oc	16.962_oc	17.062_oc	17.163_oc	17.264_oc	17.365_oc	17.466_oc	17.567_oc	17.668_oc	17.769_oc	17.87_oc	17.971_oc	18.072_oc	18.173_oc	18.274_oc	18.375_oc	18.476_oc	18.577_oc	18.678_oc	18.779_oc	18.88_oc	18.981_oc	19.082_oc	19.183_oc	19.284_oc	19.385_oc	19.486_oc	19.587_oc	19.688_oc	19.788_oc	19.889_oc	19.99_oc	20.091_oc	20.192_oc	20.293_oc	20.394_oc	20.495_oc	20.596_oc	20.697_oc	20.798_oc	20.899_oc	21.0_oc
EM10D0043278	-0.08727414570298016	-0.08680505073426538	-0.0863384770664966	-0.08587441114883614	-0.08541283950326835	-0.0849537487242083	-0.08449712547811271	-0.08404295650309278	-0.08359122860852938	-0.0831419286746898	-0.08269504365234728	-0.08225056056240179	-0.08180846649550359	-0.0813687486116782	-0.0809313941399538	-0.08049639037799046	-0.08006372469171123	-0.07963338451493564	-0.07920535734901463	-0.07877963076246781	-0.07835619239062262	-0.07793502993525521	-0.0775161311642336	-0.07709948391116236	-0.07668507607502958	-0.07627289561985541	-0.07586293057434274	-0.07545516903152961	-0.07504959914844353	-0.07464620914575766	-0.07424498730744883	-0.07384592198045739	-0.07344900157434889	-0.07305421456097763	-0.07266154947415189	-0.07227099490930115	-0.07188253952314484	-0.07149617203336316	-0.0711118812182694	-0.07072965591648427	-0.0703494850266117	-0.06997135750691669	-0.06959526237500459	-0.06922118870750242	-0.06884912563974147	-0.06847906236544218	-0.06811098813640015	-0.06774489226217419	-0.06738076410977588	-0.06701859310336099	-0.06665836872392218	-0.0663000805089839	-0.06594371805229836	-0.06558927100354353	-0.06523672906802258	-0.0648860820063651	-0.06453731963422957	-0.0641904318220079	-0.0638454084945312	-0.06350223963077731	-0.0631609152635798	-0.06282142547933857	-0.062483760417732126	-0.06214791027143114	-0.061813865285813785	-0.06148161575868249	-0.061151152039982185	-0.06082246453152028	-0.060495543686687726	-0.060170380010182115	-0.0598469640577317	-0.0595252864358214	-0.059205337801419955	-0.05888710886170864	-0.05857059037381149	-0.058255773144526915	-0.05794264803006077	-0.0576312059357608	-0.05732143781585268	-0.057013334673177284	-0.056706887558929484	-0.05640208757239828	-0.05609892586070836	-0.05579739361856313	-0.055497482087988916	-0.055199182558080796	-0.054902486364749536	-0.05460738489047009	-0.054313869564031406	-0.05402193186028754	-0.053731563299909996	-0.05344275544914172	-0.05315549991955203	-0.052869788367793165	-0.052585612495357925	-0.052302964048338826	-0.052021834817188306	-0.051742216636480466	-0.051464101384673874	-0.05118748098387577	-0.050912347399607566	-0.050638692640571435	-0.05036650875841838	-0.05009578784751738	-0.04982652204472581
EM10D0046746	-0.02097416583664084	-0.022752174901581188	-0.024670440353401078	-0.026738135372691093	-0.028964686046846076	-0.031359707109991274	-0.03393292088721107	-0.03669405715656349	-0.03965273157186402	-0.04281830028662112	-0.04619968850875734	-0.049805190924664934	-0.053642242291053334	-0.05771715703796618	-0.062034837491534475	-0.066598451344815	-0.07140908030963755	-0.0764653434936821	-0.08176300097341174	-0.08729454526354155	-0.09304879087916351	-0.09901047487504483	-0.10515988401474206	-0.11147252691077199	-0.11791887187893287	-0.1244641731122616	-0.13106840881397283	-0.1376863548252452	-0.14426781574004122	-0.15075803225221698	-0.15709827834784676	-0.1632266548808554	-0.16907907716314347	-0.17459044377249386	-0.179695962360652	-0.1843325965666337	-0.18844058711289416	-0.19196499079028143	-0.19485717432771293	-0.19707619698262502	-0.19859001672716928	-0.19937646044469925	-0.19942390847286748	-0.198731657571269	-0.19730994298164195	-0.19517961839943218	-0.1923715109222975	-0.18892548490814498	-0.18488926284359924	-0.18031706177520962	-0.17526810997225736	-0.16980511009661253	-0.16399271250289063	-0.15789605599097445	-0.1515794242471316	-0.14510505533730894	-0.13853212996854158	-0.13191595273549847	-0.1253073299801107	-0.11875213876872447	-0.11229107417174164	-0.10595955664972563	-0.0997877778696266	-0.09380086152642922	-0.0880191154741056	-0.08245835237539219	-0.07713025785298024	-0.07204278747247672	-0.0672005765519544	-0.06260534956074268	-0.058256318575960836	-0.054150562790279765	-0.0502833833312084	-0.04664862961861891	-0.043238995138740216	-0.040046281855709444	-0.03706163353644553	-0.034275739060539254	-0.03167900735799783	-0.02926171599933413	-0.027014135688926192	-0.024926633015510565	-0.02298975382141289	-0.02119428948911289	-0.01953132833055206	-0.017992294118072794	-0.01656897362947166	-0.015253534903843142	-0.014038537727569625	-0.012916937696762732	-0.011882085037684212	-0.010927719212807399	-0.010047960198746536	-0.009237297193999022	-0.008490575399402614	-0.007802981412062781	-0.007170027683616709	-0.006587536415230906	-0.006051623193747966	-0.005558680614910495	-0.005105362089621151	-0.004688565986784088	-0.004305420230519206	-0.003953267439604704	-0.0036296506721193054
```
### Examples
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

## Estimating switching time and sorting curves

### Dependencies
```
optparse
```
### Input requirements
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
### File formats
The input files should be the output values and parameters from a single chronODE run:
```
cCRE_id	10.5_oc	10.601_oc	10.702_oc	10.803_oc	10.904_oc	11.005_oc	11.106_oc	11.207_oc	11.308_oc	11.409_oc	11.51_oc	11.611_oc	11.712_oc	11.812_oc	11.913_oc	12.014_oc	12.115_oc	12.216_oc	12.317_oc	12.418_oc	12.519_oc	12.62_oc	12.721_oc	12.822_oc	12.923_oc	13.024_oc	13.125_oc	13.226_oc	13.327_oc	13.428_oc	13.529_oc	13.63_oc	13.731_oc	13.832_oc	13.933_oc	14.034_oc	14.135_oc	14.236_oc	14.337_oc	14.438_oc	14.538_oc	14.639_oc	14.74_oc	14.841_oc	14.942_oc	15.043_oc	15.144_oc	15.245_oc	15.346_oc	15.447_oc	15.548_oc	15.649_oc	15.75_oc	15.851_oc	15.952_oc	16.053_oc	16.154_oc	16.255_oc	16.356_oc	16.457_oc	16.558_oc	16.659_oc	16.76_oc	16.861_oc	16.962_oc	17.062_oc	17.163_oc	17.264_oc	17.365_oc	17.466_oc	17.567_oc	17.668_oc	17.769_oc	17.87_oc	17.971_oc	18.072_oc	18.173_oc	18.274_oc	18.375_oc	18.476_oc	18.577_oc	18.678_oc	18.779_oc	18.88_oc	18.981_oc	19.082_oc	19.183_oc	19.284_oc	19.385_oc	19.486_oc	19.587_oc	19.688_oc	19.788_oc	19.889_oc	19.99_oc	20.091_oc	20.192_oc	20.293_oc	20.394_oc	20.495_oc	20.596_oc	20.697_oc	20.798_oc	20.899_oc	21.0_oc
EM10D0043278	-0.08727414570298016	-0.08680505073426538	-0.0863384770664966	-0.08587441114883614	-0.08541283950326835	-0.0849537487242083	-0.08449712547811271	-0.08404295650309278	-0.08359122860852938	-0.0831419286746898	-0.08269504365234728	-0.08225056056240179	-0.08180846649550359	-0.0813687486116782	-0.0809313941399538	-0.08049639037799046	-0.08006372469171123	-0.07963338451493564	-0.07920535734901463	-0.07877963076246781	-0.07835619239062262	-0.07793502993525521	-0.0775161311642336	-0.07709948391116236	-0.07668507607502958	-0.07627289561985541	-0.07586293057434274	-0.07545516903152961	-0.07504959914844353	-0.07464620914575766	-0.07424498730744883	-0.07384592198045739	-0.07344900157434889	-0.07305421456097763	-0.07266154947415189	-0.07227099490930115	-0.07188253952314484	-0.07149617203336316	-0.0711118812182694	-0.07072965591648427	-0.0703494850266117	-0.06997135750691669	-0.06959526237500459	-0.06922118870750242	-0.06884912563974147	-0.06847906236544218	-0.06811098813640015	-0.06774489226217419	-0.06738076410977588	-0.06701859310336099	-0.06665836872392218	-0.0663000805089839	-0.06594371805229836	-0.06558927100354353	-0.06523672906802258	-0.0648860820063651	-0.06453731963422957	-0.0641904318220079	-0.0638454084945312	-0.06350223963077731	-0.0631609152635798	-0.06282142547933857	-0.062483760417732126	-0.06214791027143114	-0.061813865285813785	-0.06148161575868249	-0.061151152039982185	-0.06082246453152028	-0.060495543686687726	-0.060170380010182115	-0.0598469640577317	-0.0595252864358214	-0.059205337801419955	-0.05888710886170864	-0.05857059037381149	-0.058255773144526915	-0.05794264803006077	-0.0576312059357608	-0.05732143781585268	-0.057013334673177284	-0.056706887558929484	-0.05640208757239828	-0.05609892586070836	-0.05579739361856313	-0.055497482087988916	-0.055199182558080796	-0.054902486364749536	-0.05460738489047009	-0.054313869564031406	-0.05402193186028754	-0.053731563299909996	-0.05344275544914172	-0.05315549991955203	-0.052869788367793165	-0.052585612495357925	-0.052302964048338826	-0.052021834817188306	-0.051742216636480466	-0.051464101384673874	-0.05118748098387577	-0.050912347399607566	-0.050638692640571435	-0.05036650875841838	-0.05009578784751738	-0.04982652204472581
EM10D0046746	-0.02097416583664084	-0.022752174901581188	-0.024670440353401078	-0.026738135372691093	-0.028964686046846076	-0.031359707109991274	-0.03393292088721107	-0.03669405715656349	-0.03965273157186402	-0.04281830028662112	-0.04619968850875734	-0.049805190924664934	-0.053642242291053334	-0.05771715703796618	-0.062034837491534475	-0.066598451344815	-0.07140908030963755	-0.0764653434936821	-0.08176300097341174	-0.08729454526354155	-0.09304879087916351	-0.09901047487504483	-0.10515988401474206	-0.11147252691077199	-0.11791887187893287	-0.1244641731122616	-0.13106840881397283	-0.1376863548252452	-0.14426781574004122	-0.15075803225221698	-0.15709827834784676	-0.1632266548808554	-0.16907907716314347	-0.17459044377249386	-0.179695962360652	-0.1843325965666337	-0.18844058711289416	-0.19196499079028143	-0.19485717432771293	-0.19707619698262502	-0.19859001672716928	-0.19937646044469925	-0.19942390847286748	-0.198731657571269	-0.19730994298164195	-0.19517961839943218	-0.1923715109222975	-0.18892548490814498	-0.18488926284359924	-0.18031706177520962	-0.17526810997225736	-0.16980511009661253	-0.16399271250289063	-0.15789605599097445	-0.1515794242471316	-0.14510505533730894	-0.13853212996854158	-0.13191595273549847	-0.1253073299801107	-0.11875213876872447	-0.11229107417174164	-0.10595955664972563	-0.0997877778696266	-0.09380086152642922	-0.0880191154741056	-0.08245835237539219	-0.07713025785298024	-0.07204278747247672	-0.0672005765519544	-0.06260534956074268	-0.058256318575960836	-0.054150562790279765	-0.0502833833312084	-0.04664862961861891	-0.043238995138740216	-0.040046281855709444	-0.03706163353644553	-0.034275739060539254	-0.03167900735799783	-0.02926171599933413	-0.027014135688926192	-0.024926633015510565	-0.02298975382141289	-0.02119428948911289	-0.01953132833055206	-0.017992294118072794	-0.01656897362947166	-0.015253534903843142	-0.014038537727569625	-0.012916937696762732	-0.011882085037684212	-0.010927719212807399	-0.010047960198746536	-0.009237297193999022	-0.008490575399402614	-0.007802981412062781	-0.007170027683616709	-0.006587536415230906	-0.006051623193747966	-0.005558680614910495	-0.005105362089621151	-0.004688565986784088	-0.004305420230519206	-0.003953267439604704	-0.0036296506721193054
```
```
cCRE_id	        k	                    b	                  MSE	            sign_func	TYPE	    MOVE	              group	  region
EM10D0043278	-0.053382675494526265	124805.313024191	 0.02775699273405798	  1.0	upward	  0.9962693785260658	downreg	forebrain
EM10D0046746	-0.8540763957605665	  0.9343163272858892	0.011764908795083247	1.0	original	0.0	                downreg	forebrain
```
The output file is identical to the parameters input, but with columns added:
```
cCRE_id	k	b	MSE	sign_func	TYPE	MOVE	group	region	switching_time	saturation_time	minimum_time	label
EM10D0043278	-0.0533826754945263	124805.313024191	0.027756992734058	1	upward	0.996269378526066	decreasing	forebrain	-200.109828054027	-286.188688050577	490.027208851933	decelerator
EM10D0046746	-0.854076395760566	0.934316327285889	0.0117649087950832	1	original	0	decreasing	forebrain	14.6963669566587	9.31614585121227	57.8322757240436	switcher
```
### Example
```
Rscript switching.time.labels.R \
    -p example_data/oc_params.tsv \
    -v example_data/oc_vals.tsv \
    -o example_data/oc_labeled.tsv \
    -s 10.5 \
    -e 21
```

## Random Forest

### Dependencies
```
numpy
pandas
sklearn
```
### Input requirements
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
### File formats
The input file should be a tab-separated merge of the derivatives outputted by chronODE. Each row correspond to a linked gene-cCRE pair. The first two columns are expected to be indexes (i.e. the IDs of the gene-cCRE pairs). The next set of columns are the derivatives of chromatin accessibility (or other features) at each time point. The final set of columns are the derivatives of gene expression at each time point. Additions columns after the RNA columns will be ignored:
```
	cCRE_id	gene_id	10.5_oc	10.601_oc	10.702_oc	10.803_oc	10.904_oc	11.005_oc	11.106_oc	11.207_oc	11.308_oc	11.409_oc	11.51_oc	11.611_oc	11.712_oc	11.812_oc	11.913_oc	12.014_oc	12.115_oc	12.216_oc	12.317_oc	12.418_oc	12.519_oc	12.62_oc	12.721_oc	12.822_oc	12.923_oc	13.024_oc	13.125_oc	13.226_oc	13.327_oc	13.428_oc	13.529_oc	13.63_oc	13.731_oc	13.832_oc	13.933_oc	14.034_oc	14.135_oc	14.236_oc	14.337_oc	14.438_oc	14.538_oc	14.639_oc	14.74_oc	14.841_oc	14.942_oc	15.043_oc	15.144_oc	15.245_oc	15.346_oc	15.447_oc	15.548_oc	15.649_oc	15.75_oc	15.851_oc	15.952_oc	16.053_oc	16.154_oc	16.255_oc	16.356_oc	16.457_oc	16.558_oc	16.659_oc	16.76_oc	16.861_oc	16.962_oc	17.062_oc	17.163_oc	17.264_oc	17.365_oc	17.466_oc	17.567_oc	17.668_oc	17.769_oc	17.87_oc	17.971_oc	18.072_oc	18.173_oc	18.274_oc	18.375_oc	18.476_oc	18.577_oc	18.678_oc	18.779_oc	18.88_oc	18.981_oc	19.082_oc	19.183_oc	19.284_oc	19.385_oc	19.486_oc	19.587_oc	19.688_oc	19.788_oc	19.889_oc	19.99_oc	20.091_oc	20.192_oc	20.293_oc	20.394_oc	20.495_oc	20.596_oc	20.697_oc	20.798_oc	20.899_oc	21.0_oc	10.5_rna	10.601_rna	10.702_rna	10.803_rna	10.904_rna	11.005_rna	11.106_rna	11.207_rna	11.308_rna	11.409_rna	11.51_rna	11.611_rna	11.712_rna	11.812_rna	11.913_rna	12.014_rna	12.115_rna	12.216_rna	12.317_rna	12.418_rna	12.519_rna	12.62_rna	12.721_rna	12.822_rna	12.923_rna	13.024_rna	13.125_rna	13.226_rna	13.327_rna	13.428_rna	13.529_rna	13.63_rna	13.731_rna	13.832_rna	13.933_rna	14.034_rna	14.135_rna	14.236_rna	14.337_rna	14.438_rna	14.538_rna	14.639_rna	14.74_rna	14.841_rna	14.942_rna	15.043_rna	15.144_rna	15.245_rna	15.346_rna	15.447_rna	15.548_rna	15.649_rna	15.75_rna	15.851_rna	15.952_rna	16.053_rna	16.154_rna	16.255_rna	16.356_rna	16.457_rna	16.558_rna	16.659_rna	16.76_rna	16.861_rna	16.962_rna	17.062_rna	17.163_rna	17.264_rna	17.365_rna	17.466_rna	17.567_rna	17.668_rna	17.769_rna	17.87_rna	17.971_rna	18.072_rna	18.173_rna	18.274_rna	18.375_rna	18.476_rna	18.577_rna	18.678_rna	18.779_rna	18.88_rna	18.981_rna	19.082_rna	19.183_rna	19.284_rna	19.385_rna	19.486_rna	19.587_rna	19.688_rna	19.788_rna	19.889_rna	19.99_rna	20.091_rna	20.192_rna	20.293_rna	20.394_rna	20.495_rna	20.596_rna	20.697_rna	20.798_rna	20.899_rna	21.0_rna	correlation	k_rna	k_oc	region	group
0	EM10D1015674	ENSMUSG00000031138	0.0292021391502785	0.0310246038870716	0.0329472699456372	0.0349738287240955	0.037107853042689	0.0393527523719887	0.0417117227332649	0.0441876910579055	0.0467832538595565	0.0495006101588252	0.0523414887081919	0.0553070696963035	0.0583979012677808	0.0616138113782366	0.0649538157147556	0.0684160226489837	0.071997536451259	0.0756943602763145	0.079501300728531	0.0834118761198993	0.0874182308367807	0.091511058519773	0.0956795370195724	0.0999112783034576	0.104192296632836	0.108506998391806	0.112838196899277	0.117167155362931	0.121473660814789	0.125736131391604	0.129931758681014	0.134036686045857	0.138026222873155	0.141875093590069	0.145557719077199	0.149048526831373	0.152322284937744	0.155354453665117	0.158121547365217	0.160601498404518	0.162774014151848	0.164620917644638	0.166126462507241	0.167277613025171	0.168064280997973	0.168479512086641	0.168519615801813	0.168184234987419	0.16747635256352	0.166402235310595	0.164971316506698	0.163196021169346	0.16109153941233	0.158675554922904	0.155967936733888	0.152990403265893	0.149766168028246	0.146319576397016	0.142675742559718	0.138860195071462	0.134898538562459	0.130816138036583	0.126637830973107	0.122387671155973	0.118088706869063	0.113762794865966	0.10943045039261	0.105110732543756	0.10082116339114	0.096577678643045	0.0923946070836114	0.0882846756889197	0.0842590371132065	0.0803273161657045	0.0764976719370234	0.0727768723631208	0.0691703782142338	0.0656824337462482	0.0623161615352123	0.0590736593166166	0.0559560969566884	0.0529638119827631	0.0500964023856964	0.0473528156733652	0.044731433396624	0.0422301505853475	0.0398464497215207	0.0375774690388903	0.0354200650755234	0.0333708695183183	0.0314263404691099	0.0295828083327463	0.0278365165807082	0.0261836576818629	0.0246204045169811	0.0231429376078654	0.0217474684972148	0.020430259613462	0.0191876409472822	0.0180160238546252	0.0169119122861155	0.0158719117254605	0.0148927361009084	0.0139712129144637	0.013104286814022	0.110106777956121	0.110665757107439	0.111212696967868	0.111747317491891	0.112269342845492	0.112778501693311	0.113274527484222	0.113757158734829	0.114226139310393	0.114681218702713	0.115122152304437	0.11554870167934	0.115960634828052	0.116357726448758	0.116739758192374	0.117106518911731	0.11745780490427	0.117793420147805	0.118113176528875	0.118416894063248	0.11870440110814	0.118975534565726	0.119230140077524	0.119468072209283	0.11968919462597	0.119893380256522	0.120080511448008	0.120250480108893	0.1204031878411	0.120538546060592	0.120656476106234	0.120756909336691	0.120839787215165	0.120905061381795	0.12095269371356	0.120982656371557	0.120994931835567	0.120989512925815	0.120966402811892	0.120925615008818	0.120867173360249	0.120791112008868	0.120697475354025	0.120586317996709	0.120457704671977	0.120311710168975	0.120148419238732	0.119967926489908	0.119770336272723	0.119555762551316	0.119324328764784	0.119076167677209	0.11881142121697	0.118530240305679	0.118232784677089	0.117919222686343	0.117589731109949	0.117244494936893	0.116883707151288	0.116507568507014	0.116116287294772	0.115710079102021	0.11528916656625	0.114853779122072	0.114404152742603	0.113940529675624	0.113463158175013	0.11297229222793	0.112468191278268	0.111951119946846	0.111421347748845	0.110879148808986	0.110324801574927	0.109758588529364	0.109180795901321	0.1085917133771	0.10799163381134	0.107380852938659	0.106759669086309	0.106128382888285	0.105487297001309	0.104836715823091	0.104176945213281	0.103508292217472	0.102831064794643	0.102145571548377	0.10145212146221	0.100751023639413	0.100042587047531	0.099327120267949	0.0986049312507757	0.0978763270752827	0.0971416137161467	0.0964010958157099	0.0956550764624636	0.0949038569759368	0.0941477366981614	0.0933870127918608	0.0926219800454985	0.0918529306853005	0.0910801541943541	0.0903039371388652	0.0895245630016425	0.0887423120228607	0.0879574610481399	0.81216813024026	0.16939688531092	0.661612298960005	forebrain	late_upreg
2	EM10D1049003	ENSMUSG00000026085	0.284229691840175	0.281767561508326	0.278749667166499	0.275200457394488	0.271148173462744	0.266624383462629	0.261663476923738	0.256302133282552	0.250578777584522	0.244533036337194	0.238205205541685	0.231635741687187	0.22486478498282	0.217931722412557	0.210874796420489	0.203730763248036	0.196534603224406	0.189319283716667	0.182115574022282	0.174951910266073	0.167854307363119	0.160846314333842	0.153949008702044	0.147181025356532	0.140558615091995	0.134095728040167	0.127804117331557	0.121693458563799	0.115771480968807	0.110044106542912	0.104515593810409	0.0991886833124165	0.0940647423343392	0.0891439067938345	0.0844252185977047	0.0799067571337297	0.0755857638877605	0.0714587594650135	0.0675216525467517	0.0637698405299134	0.0601983017793975	0.0568016795728068	0.0535743579382671	0.0505105296804128	0.0476042569608726	0.0448495248506196	0.0422402883053237	0.0397705130341199	0.037434210739473	0.0352254692034227	0.0331384776853799	0.0311675480806464	0.02930713226844	0.0275518360547001	0.0258964300894508	0.0243358581118564	0.0228652428490516	0.0214798898679628	0.0201752896530461	0.0189471181575557	0.0177912360518061	0.0167036868691039	0.0156806942287063	0.0147186582953325	0.0138141516164785	0.012963914461993	0.0121648497750656	0.0114140178298666	0.0107086306785154	0.0100460464587394	0.0094237636234545	0.0088394151444507	0.0082907627343237	0.007775691123671	0.0072922024242853	0.0068384106035537	0.0064125360904337	0.0060129005291629	0.005637921693187	0.005286108568641	0.0049560566139742	0.0046464431999954	0.0043560232326219	0.0040836249589403	0.003828145955783	0.0035885492988494	0.0033638599094419	0.0031531610751038	0.0029555911398129	0.0027703403589057	0.0025966479135212	0.0024337990790878	0.0022811225421858	0.0021379878599995	0.002003803056525	0.0018780123496936	0.0017600940036101	0.0016495583001839	0.0015459456245351	0.0014488246586872	0.0013577906781982	0.0012724639465519	0.0011924882022946	0.0011175292340872	0.0010472735390207	1.40716563847581	1.03119604214313	0.694939955232681	0.442448607936065	0.271629967360648	0.163063994728524	0.0965787635974042	0.0567455202607257	0.0331847758865464	0.0193531240575423	0.0112684994947632	0.0065550361459686	0.0038110784518739	0.0022150490686936	0.0012871791306739	0.0007479081077319	0.0004345407798043	0.0002524626737483	0.0001466745339571	8.52132188345323e-05	4.95058080132018e-05	2.87609715061859e-05	1.67089790241285e-05	9.70723791024493e-06	5.63950626909033e-06	3.27631983124464e-06	1.90340575920652e-06	1.1057995261704e-06	6.42423441120945e-07	3.73221228023557e-07	2.16825962329006e-07	1.25966835319886e-07	7.31814739238233e-08	4.25153816458914e-08	2.46996634737227e-08	1.43494721155843e-08	8.33644489547063e-09	4.84312678630934e-09	2.81365376776536e-09	1.6346158476738101e-09	9.49644359493146e-10	5.51703350964589e-10	3.20517557823268e-10	1.86207726018624e-10	1.08177698542425e-10	6.284781690830379e-11	3.65126826343818e-11	2.1212391798533097e-11	1.2322403593215098e-11	7.15760513258276e-12	4.1603074352565895e-12	2.41741036624901e-12	1.40483165144701e-12	8.160478620092979e-13	4.73938599739813e-13	2.74981440122083e-13	1.6094501936557499e-13	9.46259661596604e-14	5.33787501413472e-14	3.31595266029582e-14	1.94104545968536e-14	1.2131534123033501e-14	5.661382590748981e-15	3.23507576614228e-15	2.4263068246067102e-15	8.087689415355689e-16	2.4263068246067102e-15	0.0	0.0	0.0	0.0	8.087689415355689e-16	8.087689415355689e-16	-1.61753788307114e-15	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0	8.087689415355689e-16	0.0	8.087689415355689e-16	0.0	8.087689415355689e-16	0.0	8.087689415355689e-16	0.0	8.087689415355689e-16	0.0	8.087689415355689e-16	8.087689415355689e-16	0.0	8.087689415355689e-16	0.0	0.0	8.087689415355689e-16	0.0	0.0	0.0	8.087689415355689e-16	0.0	0.489943751344613	5.37904223774371	0.644327248946338	forebrain	late_upreg
```
The output file will look like chronODE output, but only gene-cCRE pairs from the test set will be written to the predictions file. Note that column headers will also be replaced by generic integers:
```
combined_index	0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104
ENSMUSG00000039703-forebrain-EM10D3074493	-0.2625324998352258	-0.25478173365352913	-0.24735001435865087	-0.24022345509986845	-0.23338882387703863	-0.22683346487937184	-0.22054521936040727	-0.21451234635278596	-0.20872344393435788	-0.2031673721599803	-0.19783317905158018	-0.1927100310188509	-0.18778714855714387	-0.18305374685531725	-0.17849897899400907	-0.17411187695877983	-0.1698812834241137	-0.16579576644846877	-0.16184351163944113	-0.1580121939217417	-0.15428884500763715	-0.15065975235645024	-0.14711044691467307	-0.1436258524419062	-0.14019066802376542	-0.13679002731587842	-0.13341041941659226	-0.13004077563008787	-0.12667354676293738	-0.1233055487981739	-0.11993836797061397	-0.1165781964939453	-0.11323509665365306	-0.10992182169760253	-0.10665241271261165	-0.10344081524508278	-0.10029971923118736	-0.09723974483150323	-0.09426900711589832	-0.09139302020105688	-0.08861485898767281	-0.08593548378997674	-0.0833541418200932	-0.0808687794657065	-0.07847642203646771	-0.07617349777628665	-0.07395609807493903	-0.07182017569054806	-0.06976168825545172	-0.06777669659741505	-0.0658614276243906	-0.06401231060815621	-0.06222599429445764	-0.06049935075796555	-0.05882947052382243	-0.057213652289349104	-0.05564938962372806	-0.05413435628603767	-0.05266639125174541	-0.05124348413887118	-0.04986376144360914	-0.04852547380202888	-0.047226984365890436	-0.04596675829817516	-0.04474335334374406	-0.04355541140236351	-0.04240165101780027	-0.041280860692594505	-0.04019189293984435	-0.03913365898844064	-0.03810512406503823	-0.03710530318359633	-0.036133257380901486	-0.03518809034372858	-0.03426894537997661	-0.03337500269215228	-0.03250547691694472	-0.031659614899353414	-0.030836693673953485	-0.03003601862946141	-0.029256921835853763	-0.028498760515966563	-0.027760915645799078	-0.027042790669734574	-0.026343810318593928	-0.02566341951991922	-0.025001082391154617	-0.02435628130749799	-0.02372851603715244	-0.023117302937539032	-0.022522174206756125	-0.02194267718520398	-0.02137837370284609	-0.020828839468064717	-0.02029366349449452	-0.019772447562594563	-0.019264805713051093	-0.018770363769396302	-0.018288758887490226	-0.017819639129741716	-0.017362663062153103	-0.016917499372453463	-0.016483826507750404	-0.0160613323302765	-0.01564971378993606
ENSMUSG00000041895-midbrain-EM10D2386250	0.21748892698510805	0.16237279764122636	0.11835493752652373	0.08402078093616874	0.057404953956934464	0.03667509197989748	0.020395645660269936	0.007538533157808241	-0.0026181435389332247	-0.010604484205879774	-0.016839240428554676	-0.021671262520441673	-0.025395899659609835	-0.028260966160822854	-0.030470703602105194	-0.03219067502576487	-0.03355371001172213	-0.03466618530245821	-0.03561400042358997	-0.036467920265436306	-0.037288231188238624	-0.03812882884036442	-0.03904093173558869	-0.0400766031366697	-0.04129214112477244	-0.042751086694481986	-0.0445259617335668	-0.04669669948343552	-0.04934197976992948	-0.05251774667935265	-0.05621683868746943	-0.06030896729202473	-0.06447659609688745	-0.06818905320042179	-0.07077430710507328	-0.07161485518364824	-0.07039708582013877	-0.0672590626262324	-0.06272663264911613	-0.057486487906253145	-0.052151876985400795	-0.0471367290942729	-0.0426472168931231	-0.03873687909709961	-0.03537232202722269	-0.032483197020144274	-0.029991814287123993	-0.027827201811673676	-0.025929869936396645	-0.024251980242833975	-0.02275570987661383	-0.021411214981439147	-0.02019478294241257	-0.019087350966733606	-0.0180733885173023	-0.01714007893871171	-0.016276726235956305	-0.01547432336296724	-0.014725233479544202	-0.014022949478624034	-0.01336190792486208	-0.01273734137060457	-0.012145158391952224	-0.011581844263011393	-0.01104437751221107	-0.010530159096861695	-0.010036951884155609	-0.00956282873605077	-0.009106127890963598	-0.008665414598761535	-0.008239448148355475	-0.007827153560319449	-0.0074275973192188816	-0.007039966602711676	-0.006663551533521023	-0.006297730039814034	-0.005941954961539888	-0.005595743086235942	-0.00525866583852906	-0.004930341383659918	-0.004610427937257331	-0.004298618101757555	-0.003994634074597885	-0.003698223594959561	-0.0034091565147753118	-0.0031272218961467726	-0.0028522255516053243	-0.00258398795599447	-0.0023223424694382624	-0.0020671338200269384	-0.0018182168027893193	-0.001575455158280995	-0.0013387205999807225	-0.0011078919646357666	-0.000882854464012256	-0.0006634990201112623	-0.0004497216690704091	-0.00024142302158402192	-3.8507769976289625e-05	0.00015911576605011745	0.0003515360613230426	0.0005388387797262609	0.0007211071451451441	0.0008984222916715194	0.0010708635950900517
```
### Example
```
python randomforest.py \
  -i myfolder/rf_input.tsv \
  -p myfolder/rf_predictions.tsv \
  -s 42
 ```

## Neural Network
Coming soon!
