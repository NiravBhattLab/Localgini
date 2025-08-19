clear
initCobraToolbox(false)
load('../../InputData/TissueExpressionData.mat') %%%%%%%%%
cell_lines = expressionData.Tissue;
load('.\cons_ids.mat')
load('.\TissueInitialModel.mat')
Cmodel=removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),a)));
clear model
%% fastcore
load('LG_inputs/FASTCORE.mat')%%%%%%%%
filename = 'csm_models/Localgini/FASTCORE/';%%%%% 
for i=1:54
    core = find(core_rxn(a,i));
    m1 = fastcore(Cmodel,core,1e-8);
    model = [filename,cell_lines{i},'.mat'];
    save(fn,'model')
end

%% GIMME
load('LG_inputs/GIMME.mat')%%%%%%%%
filename = 'csm_models/Localgini/GIMME/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    model = GIMME(Cmodel,core,0,0.9);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'model')
end
   
%% iMAT
load('LG_inputs/iMAT.mat')%%%%%%%%
filename = 'csm_models/Localgini/iMAT/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    model = iMAT(Cmodel,core,1,9,1e-8,{});
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'model')
end

%% MBA
load('LG_inputs/MBA.mat')%%%%%%
filename = 'csm_models/Localgini/MBA/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    Models{i} = MBA(Cmodel,Cmodel.rxns(core==1),Cmodel.rxns(core==2),1e-8);
end
for i = 1:54
    fn = [filename,cell_lines{i},'.mat'];
    model = Models{i};
    save(fn,'model')
    clear model
end


%% mCADRE
load('LG_inputs/mCADRE.mat')%%%%%%
filename = 'csm_models/Localgini/mCADRE/';%%%%% 
for i=1:54
    ubiq = core_rxn(a,i);
    model = mCADRE(Cmodel,ubiq,zeros(numel(Cmodel.rxns),1),{},0,1/3,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'model')
    clear model
end


%% INIT
load('LG_inputs/INIT.mat')%%%%%%%%
filename = 'csm_models/Localgini/INIT/';%%%%% 
for i=1:54
    core = core_rxn(a,i);
    model = INIT(Cmodel,core,1e-8);
    fn = [filename,cell_lines{i},'.mat'];
    save(fn,'model')
end