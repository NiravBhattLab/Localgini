clear
% initCobraToolbox(false)
load('../../InputData/TissueExpressionData.mat') %%%%%%%%%

expressionDataLG.value = expressionData.valuebyTissue;
expressionDataLG.genes = expressionData.gene;
expressionDataLG.context = expressionData.Tissue;

load('.\TissueInitialModel.mat')
MeMs = {'FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE'};
ut = 75; lt = 25; ThS=1;
manual_core = {'biomass_reaction','DM_atp_c_'};
coreRxn = find(ismember(model.rxns,manual_core));

for MeM = 1:numel(MeMs)
    MeM = MeMs{MeM};
    [core_rxn,~] = GiniReactionImportance(expressionDataLG,model,MeM,ut,lt,ThS,coreRxn);
    save(['LG_inputs/',MeM],'core_rxn')
end
