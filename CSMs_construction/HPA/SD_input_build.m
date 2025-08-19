clear
load('../../InputData/TissueExpressionData.mat')
load('.\TissueInitialModel.mat')
modelData= getModelData(expressionData,model);
spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);

biomass_id = find(strcmp(model.rxns,'biomass_reaction'));
atp_demand_id = find(strcmp(model.rxns,'DM_atp_c_'));
coreRxn=[biomass_id atp_demand_id];


enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3 -2 -1 0 1 2 2.5 3 4 5]; % bins  
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,19,distMethod,linkageMethod);

% FASTCORE
core = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]);
core(coreRxn,:) = 1;
save('SD_inputs\FASTCORE','core')
clear core

% GIMME
core = getGIMMEscores(clustObj,edgeX,model);
core(core==-inf)=min(min(core(core~=-inf)));
core(coreRxn,:)=1;
core(find(sum(model.rxnGeneMat,2)==0),:)=0;
save('SD_inputs\GIMME','core')
clear core

% MBA
[H,M] = getMBAsets(clustObj,edgeX,model,0.1);
H(coreRxn,:)=1;
save('SD_inputs\MBA','H','M')

% INIT
core = getINITweights(clustObj,edgeX,model);
core(core==-inf)=min(min(core(core~=-inf)));
core(coreRxn,:)=1;
core(find(sum(model.rxnGeneMat,2)==0),:)=0;
save('SD_inputs\INIT','core')

% mCADRE
[core,~] = getUbiquityScore(clustObj,edgeX,model);
core(coreRxn,:)=1;
core(find(sum(model.rxnGeneMat,2)==0),:)=-1;
save('SD_inputs\mCADRE','core')

% iMAT
core = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]);
core(coreRxn,:)=1;
core = double(core);
core(find(sum(model.rxnGeneMat,2)==0),:)=0.5;
save('SD_inputs\iMAT','core')
