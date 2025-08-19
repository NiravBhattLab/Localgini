# Details of the folders and code files

1) LG_input_build.m: Code file to get the inputs required for distinct MeMs from Localgini thresholding method
2) LT2_input_build.m: Code file to get the inputs required for distinct MeMs from LocalT2 thresholding method
3) SD_input_build.m: Code file to get the inputs required for distinct MeMs from StanDep thresholding method
4) LG_model_build.m: Code file to build the context-specific models using Localgini thresholding method
5) LT2_model_build.m: Code file to build the context-specific models using LocalT2 thresholding method
6) SD_model_build.m: Code file to build the context-specific models using StanDep thresholding method
7) TissueInitialModel.mat: Recon2.2 GEM with gene IDs that are similar to the gene expression data
8) cons_ids.mat: Indices of reactions that are consistent in TissueInitialModel.mat
9) csm_models: the constructed context-specific models are stored in this folder
10) LG_inputs: The .mat files that are obtained after running LG_input_build.m are stored here
11) LT2_inputs: The .mat files that are obtained after running LT2_input_build.m are stored here
12) SD_inputs: The .mat files that are obtained after running SD_input_build.m are stored here