clear
load('../../InputData/TissueExpressionData.mat')
load('.\TissueInitialModel.mat')
modelData= getModelData(expressionData,model);

biomass_id = find(strcmp(model.rxns,'biomass_reaction'));
atp_demand_id = find(strcmp(model.rxns,'DM_atp_c_'));
coreRxn=[biomass_id atp_demand_id];

% FASTCORE
rxnTisMat = getFCcore_lt2(modelData,model,25,75);
rxnTisMat(coreRxn,:) =1;
core = rxnTisMat;
save('LT2_inputs\FASTCORE','core')

% GIMME
gimme_scores = getGIMMEscores_lt2(modelData,model,25,75);
gimme_scores(coreRxn,:) =10*log(2);
core = gimme_scores;
ids =find(sum(model.rxnGeneMat,2)==0);
core(ids,:)=0;
save('LT2_inputs\GIMME','core')

% MBA
[H,M] = getMBAsets_lt2(modelData,model,25,75);
H(coreRxn,:)=1;
save('LT2_inputs\MBA','H','M')

% INIT
weights = getINITweights_lt2(modelData,model,25,75);
weights(coreRxn,:) =max(max(weights));
core = weights;
ids =find(sum(model.rxnGeneMat,2)==0);
core(ids,:)=0;
save('LT2_inputs\INIT','core')
