clear
initCobraToolbox(false)
load('../../InputData/CancerExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

load('generic_models\UACC257.mat') %%%%%%%%%%
modelData = getModelData(expressionData,model);

spec = getSpecialistEnzymes(model);  
prom = getPromEnzymes(model);
enzymeData = comparePromiscuousSpecific(spec,prom,modelData);
edgeX = [-3 -2 -1 0 1 2 2.5 3 4]; % bins  
distMethod = 'euclidean'; % distance method  
linkageMethod = 'complete'; % linkage metric for hierarchical clustering
clustObj = geneExprDist_hierarchy(enzymeData,[],edgeX,26,distMethod,linkageMethod);
core_rxn = models4mClusters1(clustObj,enzymeData.Tissue,model,edgeX,[],[],false,0,[1 1]);

% fastcore and iMAT
save('SD_inputs/inFASTCORE','core_rxn')
save('SD_inputs/iniMAT','core_rxn')
% GIMME
core_rxn = getGIMMEscores(clustObj,edgeX,model);
core_rxn(core_rxn==-inf)=min(min(core_rxn(core_rxn~=-inf)));
save('SD_inputs/inGIMME','core_rxn')
% mCADRE
[core_rxn,~] = getUbiquityScore(clustObj,edgeX,model);
save('SD_inputs/inmCADRE','core_rxn')
% MBA
[H,M] = getMBAsets(clustObj,edgeX,model,0.1);
save('SD_inputs/inMBA','H','M')
% INIT
core_rxn = getINITweights(clustObj,edgeX,model);
core_rxn(core_rxn==-inf)=min(min(core_rxn(core_rxn~=-inf)));
save('SD_inputs/inINIT','core_rxn')

