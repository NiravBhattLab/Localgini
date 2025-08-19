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
manual_core = {'biomass_reaction','DM_atp_c_'};

%% fastcore
load('LT2_inputs/inFASTCORE.mat')%%%%%%%%
filename = 'csm_models/LocalT2/FASTCORE/';%%%%% 
core_rxn(ismember(gen_mod.rxns,manual_core),:)=1;
for i=1:44
   load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = find(core_rxn(ismember(gen_mod.rxns,model.rxns),i));
    m1 = fastcore(model,core,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end

%% GIMME
load('LT2_inputs/inGIMME.mat')%%%%%%%%
filename = 'csm_models/LocalT2/GIMME/';%%%%% 
core_rxn(ismember(gen_mod.rxns,manual_core),:)=10*log(2);
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=5*log(2);%%%%%%%%%%%%%%%%%%%%%%
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = GIMME(model,core,5*log(2));
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end


%% iMAT
load('LT2_inputs/iniMAT.mat')%%%%%%%%
filename = 'csm_models/LocalT2/iMAT/';%%%%% 
core_rxn(ismember(gen_mod.rxns,manual_core),:)=10*log(2);
% core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=0.5;
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = iMAT(model,core,5*log(2),5*log(2),1e-8,{});
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end

%% MBA
load('LT2_inputs/inMBA.mat')%%%%%%
filename = 'csm_models/LocalT2/MBA/';%%%%% 
H(ismember(gen_mod.rxns,manual_core),:)=1;
parfor i=1:44
    initCobraToolbox(false)
    model = load(['consistent_genericmodels/',cell_lines{i},'.mat']);%%%%%%%
    model = model.model;
    curr_H = model.rxns(find(H(ismember(gen_mod.rxns,model.rxns),i)));
    curr_M = model.rxns(find(M(ismember(gen_mod.rxns,model.rxns),i)));
    Models{i} = MBA(model,curr_M,curr_H,1e-8);
end
for i = 1:44
    fn = [filename,cell_lines{i},'.mat'];
    m1 = Models{i};
    save(fn,'m1')
end
%% mCADRE
load('LT2_inputs/inmCADRE.mat')%%%%%%
filename = 'csm_models/LocalT2/mCADRE/';%%%%% 
core_rxn(ismember(gen_mod.rxns,manual_core),:)=10*log(2);
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=-1; %%%%%%%%%%%%%%
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
load('LT2_inputs/inINIT.mat')%%%%%%%%
filename = 'csm_models/LocalT2/INIT/';%%%%% 
core_rxn(ismember(gen_mod.rxns,manual_core),:)=max(max(core_rxn));
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=0;%%%%%%%%%%%%%%%%
for i=1:44
    load(['consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'m1')
end



