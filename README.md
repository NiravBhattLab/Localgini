# Localgini
Codes for "Localgini: A gini coefficient based thresholding algorithm to build constraint based
context-specific metabolic models"   

Authors: Pavan Kumar S and Nirav P Bhatt 

### Requirements
1. MATLAB
2. [COBRA Toolbox](http://opencobra.github.io/cobratoolbox/)

### Note
1. Gene essentiality predictions uses ```fastSLgenes.m``` from https://github.com/RamanLab/FastSL

### To get a context specific model using Localgini follow the steps below
1. mRNA expression data has to be available as a matlab structure with fields:   
  	.value : mRNA gene expression matrix with dimension N_genes*N_samples <br>
	.genes : cell array with geneIDs in the same format as model.genes <br>
	.context : cell array with names of the samples <br>

2. genome-scale model has to be available as a COBRA model structure.


**initializing COBRA toolbox**  
```
initCobraToolbox
```

**Loading the genome scale model**  
```
model = readCbModel('.\Data\ConsRecon.mat');
```

**Loading gene expression data**  
```
load('.\Data\geneExpression.mat');
```

**Choice of model extraction method (has to be anyone of ['FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE'])**  
```
MeM = 'FASTCORE';
```

**Specifying contexts for which models have to be built (has to be in the same format as geneExpression.context)**  
```
contexts = {'C1','C2','C4','C20','C35'};
```

**Specifying upper threshold percentile (Global threshold above which all the reactions are considered active)**  
```
ut = 75;
```

**Specifying lower threshold percentile (Global threshold below which all the reactions are considered inactive)**  
```
lt = 25;
```

**Specifying where threshold has to be implied (1==> implied to genes; 2==> implied to enzymes; 3==> implied to reactions)**  
```
ThS = 1; % impliying at gene level
```

**Reactions that are manually given higher importance**
```
biomass_id = find(strcmp(model.rxns,'biomass_reaction'));
atp_demand_id = find(strcmp(model.rxns,'DM_atp_c_'));
coreRxn=[biomass_id atp_demand_id];
```

**Tolerance level above which reactions are considered as expressed**
```
tol = 1e-4;
```

**Folder path to save the built models**
```
filename = './';
```

***Reaction ids of consistent reactions in the model***
```
cons_mod_rxn_id = [1:numel(model.rxns)];
```

***Extracting context-specific models***
```
[Models,RxnImp] = buildContextmodels(geneExpression,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,cons_mod_rxn_id,tol);
```

__________________________________________________________________________

### Description of other .m functions

```GeneEssentialityAcc.m``` : To compute accuracy, pvalue and essential genes for the given context and CRISPR-CAS9 scores

```SelfConsistency.m``` : To compute Jaccard similarity betweeen core reaction list and model's reactions and fractional contibution from MeM

```GiniReactionImportance.m``` : Returns MeM specific inputs required to build context-specific models

__________________________________________________________________________

### Acknowledgement
* [pCoE for Network Systems Learning, Control and Evolution](https://ioe.iitm.ac.in/project/network-systems/)
* [Centre for Integrative Biology and Systems medicinE](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)
