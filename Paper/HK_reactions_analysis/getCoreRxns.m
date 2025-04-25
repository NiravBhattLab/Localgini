clear
% loading the models
recon22 = load('./Recon2_2'); % recon2.2 model
recon22_cancer = recon22.model; % recon2.2 to for cancer
recon22_tissue = getTissuerecon22(recon22.model); % recon2.2 to for tissue

recon3d = load('./Recon3D'); % recon3d model
recon3d_cancer = recon3d.model; 
recon3d_cancer.genes = regexprep(recon3d_cancer.genes, '_\w*', ''); % recon3d for cancer
recon3d_tissue = getTissueRecon3d(recon3d.model);

ihuman = load('./iHuman'); % iHuman model
ihuman_cancer = getCancerIhuman(ihuman.model); % ihuman for cancer
ihuman_tissue = ihuman.model; % ihuman for tissue
%% For cancer expression data
load('./CancerExpressionData')

% localgini thresholding
addpath('../Localgini')
MeM = 'FASTCORE';
ut = 75;
lt = 25;
ThS = 1;
coreRxn={};

expressionDataLG.value = expressionData.valuebyTissue;
expressionDataLG.genes = expressionData.gene;
expressionDataLG.context = expressionData.Tissue;
[LG_core_recon22,~] = GiniReactionImportance(expressionDataLG,recon22_cancer,MeM,ut,lt,ThS,coreRxn);
[LG_core_recon3d,~] = GiniReactionImportance(expressionDataLG,recon3d_cancer,MeM,ut,lt,ThS,coreRxn);
[LG_core_ihuman,~] = GiniReactionImportance(expressionDataLG,ihuman_cancer,MeM,ut,lt,ThS,coreRxn);
rmpath('../Localgini')

% localt2 thresholding
addpath('../LocalT2')
modelData = getModelData(expressionData,recon22_cancer);
LT2_core_recon22 = getFCcore_lt2(modelData,recon22_cancer,25,75);
modelData = getModelData(expressionData,recon3d_cancer);
LT2_core_recon3d = getFCcore_lt2(modelData,recon3d_cancer,25,75);
modelData = getModelData(expressionData,ihuman_cancer);
LT2_core_ihuman = getFCcore_lt2(modelData,ihuman_cancer,25,75);
rmpath('../LocalT2')

% Standep thresholding
addpath('../StanDep')

modelData= getModelData(expressionData,recon22_cancer);
spec = getSpecialistEnzymes(recon22_cancer);  
prom = getPromEnzymes(recon22_cancer);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,26,distMethod,linkageMethod);
SD_core_recon22 = models4mClusters1(clustObj,enzymeData.Tissue,recon22_cancer,edgeX,[],[],false,0,[1 1]);

modelData= getModelData(expressionData,recon3d_cancer);
spec = getSpecialistEnzymes(recon3d_cancer);  
prom = getPromEnzymes(recon3d_cancer);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,26,distMethod,linkageMethod);
SD_core_recon3d = models4mClusters1(clustObj,enzymeData.Tissue,recon3d_cancer,edgeX,[],[],false,0,[1 1]);

modelData= getModelData(expressionData,ihuman_cancer);
spec = getSpecialistEnzymes(ihuman_cancer);  
prom = getPromEnzymes(ihuman_cancer);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,26,distMethod,linkageMethod);
SD_core_ihuman = models4mClusters1(clustObj,enzymeData.Tissue,ihuman_cancer,edgeX,[],[],false,0,[1 1]);
rmpath('../StanDep')

save('core_rxns_cancer_data','LG_core_recon22','LG_core_recon3d','LG_core_ihuman',...
    'LT2_core_recon22','LT2_core_recon3d','LT2_core_ihuman',...
    'SD_core_recon22','SD_core_recon3d','SD_core_ihuman')
clearvars -except recon22 recon22_cancer recon22_tissue recon3d recon3d_cancer recon3d_tissue ihuman ihuman_cancer ihuman_tissue

%% For HPA data
load('./TissueExpressionData')

% localgini thresholding
addpath('../Localgini')
MeM = 'FASTCORE';
ut = 75;
lt = 25;
ThS = 1;
coreRxn={};

expressionDataLG.value = expressionData.valuebyTissue;
expressionDataLG.genes = expressionData.gene;
expressionDataLG.context = expressionData.Tissue;
[LG_core_recon22,~] = GiniReactionImportance(expressionDataLG,recon22_tissue,MeM,ut,lt,ThS,coreRxn);
[LG_core_recon3d,~] = GiniReactionImportance(expressionDataLG,recon3d_tissue,MeM,ut,lt,ThS,coreRxn);
[LG_core_ihuman,~] = GiniReactionImportance(expressionDataLG,ihuman_tissue,MeM,ut,lt,ThS,coreRxn);
rmpath('../Localgini')

% localt2 thresholding
addpath('../LocalT2')
modelData = getModelData(expressionData,recon22_tissue);
LT2_core_recon22 = getFCcore_lt2(modelData,recon22_tissue,25,75);
modelData = getModelData(expressionData,recon3d_tissue);
LT2_core_recon3d = getFCcore_lt2(modelData,recon3d_tissue,25,75);
modelData = getModelData(expressionData,ihuman_tissue);
LT2_core_ihuman = getFCcore_lt2(modelData,ihuman_tissue,25,75);
rmpath('../LocalT2')

% Standep thresholding
addpath('../StanDep')
modelData= getModelData(expressionData,recon22_tissue);
spec = getSpecialistEnzymes(recon22_tissue);  
prom = getPromEnzymes(recon22_tissue);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4, 5];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,19,distMethod,linkageMethod);
SD_core_recon22 = models4mClusters1(clustObj,enzymeData.Tissue,recon22_tissue,edgeX,[],[],false,0,[1 1]);

modelData= getModelData(expressionData,recon3d_tissue);
spec = getSpecialistEnzymes(recon3d_tissue);  
prom = getPromEnzymes(recon3d_tissue);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4, 5];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,19,distMethod,linkageMethod);
SD_core_recon3d = models4mClusters1(clustObj,enzymeData.Tissue,recon3d_tissue,edgeX,[],[],false,0,[1 1]);

modelData= getModelData(expressionData,ihuman_tissue);
spec = getSpecialistEnzymes(ihuman_tissue);  
prom = getPromEnzymes(ihuman_tissue);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3, -2, -1, 0, 1, 2, 2.5, 3, 4, 5];
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,19,distMethod,linkageMethod);
SD_core_ihuman = models4mClusters1(clustObj,enzymeData.Tissue,ihuman.model,edgeX,[],[],false,0,[1 1]);
rmpath('../StanDep')

save('core_rxns_tissue_data','LG_core_recon22','LG_core_recon3d','LG_core_ihuman',...
    'LT2_core_recon22','LT2_core_recon3d','LT2_core_ihuman',...
    'SD_core_recon22','SD_core_recon3d','SD_core_ihuman')

function model = getCancerIhuman(model)
    filename = './gecko_id_to_entrez.txt';% gene symbol to entrz id
    opts = detectImportOptions(filename);
    symtoentrez = readtable(filename,opts);
    temp = {};
    for i=1:numel(model.geneShortNames)
        g = symtoentrez.Var2(ismember(symtoentrez.Var1,model.geneShortNames{i}));
        temp{i,1} = g;
    end
    model.genes = cellfun(@(x)num2str(x),temp,'UniformOutput',false);
end
function model = getTissuerecon22(model)
    filename = './entrez_to_ensembl.txt';% gene symbol to entrz id
    opts = detectImportOptions(filename);
    ent2ens = readtable(filename,opts);
    ent2ens.From = arrayfun(@(x)num2str(x),ent2ens.From,'UniformOutput',false);
    temp = {};
    for i= 1:numel(model.genes)
        g = ent2ens.To(ismember(ent2ens.From,model.genes{i}));
        if ~isempty(g)
            temp{i,1} = g;
        else
            temp{i,1} = {''};
        end
    end
    model.genes = cellfun(@(x)x{1},temp,'UniformOutput',false);
end

function model = getTissueRecon3d(model)
    filename = './ensembl_to_symbol.txt';
    opts = detectImportOptions(filename);
    ens2sym = readtable(filename,opts);
    temp = {};
    for i=1:numel(model.geneisrefseq_nameID)
        g = ens2sym.ENSEMBL(ismember(ens2sym.SYMBOL,model.geneisrefseq_nameID{i}));
        if ~isempty(g)
            temp{i,1} = g;
        else
            temp{i,1} = {''};
        end
    end
    model.genes = cellfun(@(x)x{1},temp,'UniformOutput',false);
end