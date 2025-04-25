clear

%% For cancer data
load('core_rxns_cancer_data.mat')
% For Recon2.2
modelName = 'Recon2_2';
nHK_LG_recon22 = getHKinCore(LG_core_recon22,modelName);
nHK_LT2_recon22 = getHKinCore(LT2_core_recon22,modelName);
nHK_SD_recon22 = getHKinCore(SD_core_recon22,modelName);

% For Recon3D
modelName = 'Recon3D';
nHK_LG_recon3d = getHKinCore(LG_core_recon3d,modelName);
nHK_LT2_recon3d = getHKinCore(LT2_core_recon3d,modelName);
nHK_SD_recon3d = getHKinCore(SD_core_recon3d,modelName);

% For iHuman
modelName = 'iHuman';
nHK_LG_ihuman = getHKinCore(LG_core_ihuman,modelName);
nHK_LT2_ihuman = getHKinCore(LT2_core_ihuman,modelName);
nHK_SD_ihuman = getHKinCore(SD_core_ihuman,modelName);


load('CancerExpressionData.mat')
contexts = expressionData.Tissue;

save('HK_in_core_cancer','nHK_LG_recon22','nHK_LT2_recon22','nHK_SD_recon22',...
    'nHK_LG_recon3d','nHK_LT2_recon3d','nHK_SD_recon3d',...
    'nHK_LG_ihuman','nHK_LT2_ihuman','nHK_SD_ihuman','contexts')


%% Tissue data
clear
load('core_rxns_tissue_data.mat')
% For Recon2.2
modelName = 'Recon2_2';
nHK_LG_recon22 = getHKinCore(LG_core_recon22,modelName);
nHK_LT2_recon22 = getHKinCore(LT2_core_recon22,modelName);
nHK_SD_recon22 = getHKinCore(SD_core_recon22,modelName);

% For Recon3D
modelName = 'Recon3D';
nHK_LG_recon3d = getHKinCore(LG_core_recon3d,modelName);
nHK_LT2_recon3d = getHKinCore(LT2_core_recon3d,modelName);
nHK_SD_recon3d = getHKinCore(SD_core_recon3d,modelName);

% For iHuman
modelName = 'iHuman';
nHK_LG_ihuman = getHKinCore(LG_core_ihuman,modelName);
nHK_LT2_ihuman = getHKinCore(LT2_core_ihuman,modelName);
nHK_SD_ihuman = getHKinCore(SD_core_ihuman,modelName);
load('TissueExpressionData.mat')
contexts = expressionData.Tissue;
save('HK_in_core_tissue','nHK_LG_recon22','nHK_LT2_recon22','nHK_SD_recon22',...
    'nHK_LG_recon3d','nHK_LT2_recon3d','nHK_SD_recon3d',...
    'nHK_LG_ihuman','nHK_LT2_ihuman','nHK_SD_ihuman','contexts')


function nHK = getHKinCore(coreRxns,modelName)
    load(['HK_rxns_',modelName])
    load(modelName)
    ids = ismember(model.rxns,HK_rxns);
    nHK = sum(coreRxns(ids,:),1);
end