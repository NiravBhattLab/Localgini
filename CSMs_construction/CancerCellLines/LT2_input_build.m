clear
initCobraToolbox(false)
load('../../InputData/CancerExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

load('generic_models\UACC257.mat') %%%%%%%%%%
modelData = getModelData(expressionData,model);

lowerThs=25;
upperThs=75;
thrRxnData = getLocalT2_case(modelData,model,lowerThs,upperThs);

% fastcore
core_rxn = thrRxnData.value> 5*log(2);
save('LT2_inputs/inFASTCORE','core_rxn')

% GIMME and iMAT and mCADRE
core_rxn = thrRxnData.value;
save('LT2_inputs/inGIMME','core_rxn')
save('LT2_inputs/iniMAT','core_rxn')
save('LT2_inputs/inmCADRE','core_rxn')

% INIT
core_rxn = thrRxnData.value;
core_rxn(core_rxn<5*log(2)) = -8;
core_rxn(core_rxn~=-8) = core_rxn(core_rxn~=-8)/(5*log(2));
save('LT2_inputs/inINIT','core_rxn')

% MBA 
core_rxn = thrRxnData.value;
up_thr = prctile(core_rxn(core_rxn~=-2),75,'all');
lw_thr = 5*log(2);
H = core_rxn>=up_thr;
M = core_rxn>=lw_thr & core_rxn<up_thr;
save('LT2_inputs/inMBA','H','M')
