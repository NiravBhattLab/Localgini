clear
% initCobraToolbox(false)
load('../../InputData/CancerExpressionData.mat') %%%%%%%%%

expressionDataLG.value = expressionData.valuebyTissue;
expressionDataLG.genes = expressionData.gene;
expressionDataLG.context = expressionData.Tissue;

cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

load('generic_models\UACC257.mat') %%%%%%%%%%
MeMs = {'FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE'};
ut = 75; lt = 25; ThS=1;
manual_core = {'biomass_reaction','DM_atp_c_'};
coreRxn = find(ismember(model.rxns,manual_core));

for MeM = 1:numel(MeMs)
    MeM = MeMs{MeM};
    [core_rxn,~] = GiniReactionImportance(expressionDataLG,model,MeM,ut,lt,ThS,coreRxn);
    save(['LG_inputs/in',MeM],'core_rxn')
end
