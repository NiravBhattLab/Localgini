# Localgini
Codes for "Modelling Reliable Metabolic Phenotypes by Analysing The Context-Specific Transcriptomics Data"   

Authors: Pavan Kumar S and Nirav P Bhatt 

# To generate a context-specific metabolic model using Localgini

### Requirements
1. MATLAB
2. [COBRA Toolbox](http://opencobra.github.io/cobratoolbox/)


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
model = readCbModel('GEM.mat');
```

**Loading gene expression data**  
```
load('geneExpression.mat');
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
# To reproduce the results generated 

### Requirements

1) Download the human GEMs [iHuman](https://github.com/SysBioChalmers/Human-GEM/tree/main/model), [Recon2.2](https://www.ebi.ac.uk/biomodels/MODEL1603150001#Overview) and [Recon3D](http://bigg.ucsd.edu/models/Recon3D) in .mat formats and store them at the folder InputData. These models has to be flux consistent models.

2) Download the NCI60 and HPA gene-expression data and store them at the folder InputData. The data must be in the format as defined above.

### Details of the folders

1) HK_reaction_analysis
 - To get the housekeeping reactions from the housekeeping genes for all the three GEMs.
 - To get the core reactions for cancer cell line data and HPA data
 - To get the count of housekeeping reactions rectified by each of the thresholding methods
2) CSMs_construction
 - To construct the context-spcific metabolic models using three distict thresholding methods and six different MeMs for the two gene expression datasets
3) HK_reactions_in_models
 - To get the fraction of housekeeping reactions rectified in each of the context-specific metabolic models built.
4) Cancer_hallmark_gene_recovery
 - To get the fraction of hall mark genes captured by the context-specific models
5) Pathway_enrichment
 - To do enrichment analysis on the core reactions derived from the three GEMs model to identify the enriched pathways
6) SelfConsistency_analysis
 - To do the self consistency analysis on top of the CSMs to get the fractional contribution of reactions added by the MeMs in both the gene expression data. Further to perform hypothesis test to get the pvalue on how different the distribution of fractional contribution of localgini derived models compared to the others
7) variance_in_models
 - To generate metabolic tasks report on the models built using three different thresholding methods
 - To do PCA on reaction content matrix and metabolic task matrix to get the variance contribution of each of the factor in the final CSM model
8) Visualisations
 -  All the results generated are stored in this folder as .mat, .xlsx, or .csv formats
 -  Python codes to get all the plots presented in the main text and supplementary text

### Acknowledgement
* [Centre for Integrative Biology and Systems medicinE](https://ibse.iitm.ac.in/)
* [Robert Bosch Centre for Data Science and Artificial Intelligence (RBCDSAI)](https://rbcdsai.iitm.ac.in/)
