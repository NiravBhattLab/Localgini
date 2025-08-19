clear
initCobraToolbox(false)
load('../../InputData/CancerExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

load('./generic_models/786O.mat')%%%%%%
gen_mod = model;
clear model

%% fastcore
load('LG_inputs/inFASTCORE.mat')%%%%%%%%
filename = 'csm_models/Localgini/FASTCORE/';%%%%% 
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = find(core_rxn(ismember(gen_mod.rxns,model.rxns),i));
    m1 = fastcore(model,core,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end

%% GIMME
load('LG_inputs/inGIMME.mat')%%%%%%%%
filename = 'csm_models/Localgini/GIMME/';%%%%% 
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = GIMME(model,core,0,0.9);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end
   
%% iMAT
load('LG_inputs/iniMAT.mat')%%%%%%%%
filename = 'csm_models/Localgini/iMAT/';%%%%% 
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = iMAT(model,core,1,9,1e-8,{});
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end

%% MBA
load('LG_inputs/inMBA.mat')%%%%%%
filename = 'csm_models/Localgini/MBA/';%%%%% 
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat']);
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    Models{i} = MBA(model,model.rxns(core==1),model.rxns(core==2),1e-8);
end
for i = 1:44
    fn = [filename,cell_lines{i},'.mat'];
    m1 = Models{i};
    save(fn,'m1')
end


%% mCADRE
load('LG_inputs/inmCADRE.mat')%%%%%%
filename = 'csm_models/Localgini/mCADRE/';%%%%% 
for i=1:44
    model = load(['consistent_genericmodels/',cell_lines{i},'.mat']);%%%%%%%
    model = model.model;
    ubiq = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = mCADRE(model,ubiq,zeros(numel(model.rxns),1),{},0,1/3,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
    clear m1
end


%% INIT
load('LG_inputs/inINIT.mat')%%%%%%%%
filename = 'csm_models/Localgini/INIT/';%%%%% 
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end