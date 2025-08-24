clear
% for tissue data
load('../InputData/TissueExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;

load('../CSMs_construction/HPA/TissueInitialModel.mat')
load('../CSMs_construction/HPA/cons_ids.mat')
gen_mod = model;
clear model
load('lg_tissue_rxn_imp.mat') %%%%%%%%
manual_core = {'biomass_reaction','DM_atp_c_'};
coreRxn = find(ismember(gen_mod.rxns,manual_core));
Cmodel=removeRxns(gen_mod,gen_mod.rxns(setdiff(1:numel(gen_mod.rxns),a)));


% for cutoff = -2
co=2;%%%%% 
eps=0.001;
core_rxn = RxnImp;
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=-eps;core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');
filename = 'GIMME_tissue_cell_line/co_2/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    model = GIMME(Cmodel,core,0,0.9);
    save(fn,'model')
end



% for cutoff = -6
co=6;
eps=0.001;
core_rxn = RxnImp;
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=-eps;core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');
filename = 'GIMME_tissue_cell_line/co_6/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    model = GIMME(Cmodel,core,0,0.9);
    save(fn,'model')
end



% for cutoff = -14
co=14;
eps=0.001;
core_rxn = RxnImp;
core_rxn(find(sum(gen_mod.rxnGeneMat,2)==0),:)=-eps;core_rxn(core_rxn<-co)=-co;
core_rxn(coreRxn,:)=max(core_rxn,[],'all');
filename = 'GIMME_tissue_cell_line/co_14/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    fn = [filename,cell_lines{i},'.mat'];
    if isfile(fn)
        continue; % Skip this iteration
    end
    model = GIMME(Cmodel,core,0,0.9);
    save(fn,'model')
end