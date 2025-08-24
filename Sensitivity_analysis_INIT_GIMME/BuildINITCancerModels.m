clear
% for cancer cell line data
load('../InputData/CancerExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');%%%%%%%%%
cell_lines{strcmp(cell_lines,'x786O')}='786O';%%%%%%%%%
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';%%%%%%%%%

load('../CSMs_construction/CancerCellLines/generic_models/786O.mat')
gen_mod = model;
clear model
load('lg_cancer_rxn_imp.mat') %%%%%%%%
manual_core = {'biomass_reaction','DM_atp_c_'};
coreRxn = find(ismember(gen_mod.rxns,manual_core));

% for cutoff = -2
co=2;%%%%% 
core_rxn = RxnImp;
core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');

filename = 'INIT_cancer_cell_line/co_2/';%%%%% 
for i=1:44
    load(['../CSMs_construction/CancerCellLines/consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    save(fn,'m1')
end


% for cutoff = -6
co=6;%%%%% 
core_rxn = RxnImp;
core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');

filename = 'INIT_cancer_cell_line/co_6/';%%%%% 
for i=1:44
    load(['../CSMs_construction/CancerCellLines/consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    save(fn,'m1')
end


% for cutoff = -10
co=10;%%%%% 
core_rxn = RxnImp;
core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');

filename = 'INIT_cancer_cell_line/co_10/';%%%%% 
for i=1:44
    load(['../CSMs_construction/CancerCellLines/consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    save(fn,'m1')
end


% for cutoff = -14
co=14;%%%%% 
core_rxn = RxnImp;
core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');

filename = 'INIT_cancer_cell_line/co_14/';%%%%% 
for i=1:44
    load(['../CSMs_construction/CancerCellLines/consistent_genericmodels/',cell_lines{i},'.mat'])%%%%%%%
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    core = core_rxn(ismember(gen_mod.rxns,model.rxns),i);
    m1 = INIT(model,core,1e-8);
    save(fn,'m1')
end