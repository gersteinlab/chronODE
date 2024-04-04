# *chronODE*: A framework to integrate time-series multi-omics data based on ordinary differential equations combined with machine learning
Beatrice Borsari, Mor Frank, Eve S. Wattenberg, Susanna X. Liu, Xuezhu Yu, Mark Gerstein  
[bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2023.12.13.571513v1)  

An ODE-based  pipeline for interpolating values and derivatives from time-series data, and machine learning models that predicted gene expression from linked chromatin features.

![](https://github.com/gersteinlab/chronODE/blob/main/figure1.png)

***

## *chronODE*

### Requirements

```Libraries should go here?```
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
The parameters output will be tab-separated and have a row for each element and have columns for the k and b parameters, mean squared error (MSE), ????, the vertical move applied to the data, user-specified group, and user-specified region:  
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

## Switching time classifier

### Requirements

```Libraries should go here?```
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
cCRE_id	k	b	MSE	sign_func	TYPE	MOVE	group	region	cekt_start	cekt_end	ratio_start	ratio_end	switching_time	saturation_time	minimum_time	label
EM10D0043278	-0.0533826754945263	124805.313024191	0.027756992734058	1	upward	0.996269378526066	decreasing	forebrain	1.63492053276479	0.933397715159376	1.30997670944337e-05	7.47882996758685e-06	-200.109828054027	-286.188688050577	490.027208851933	decelerator
EM10D0046746	-0.854076395760566	0.934316327285889	0.0117649087950832	1	original	0	decreasing	forebrain	33.652176943635	0.00428890249360284	36.017969461576	0.00459041800763746	14.6963669566587	9.31614585121227	57.8322757240436	switcher
```

## Random Forest
### Requirements
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
The input file should be a tab-separated merge of the derivatives outputted by chronODE. Each row correspond to a linked gene-cCRE pair. The first two columns are expected to be indexes (i.e. the IDs of the gene-cCRE pairs). The next set of columns are the derivatives of chromatin accessibility (or other features) at each time point. The final set of columns are the derivatives of gene expression at each time point:
```
tk
```
The output file will look like chronODE output, but only gene-cCRE pairs from the test set will be written to the predictions file:
```
tk
```
## Neural Net
