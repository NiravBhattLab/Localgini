clear
load('../../InputData/TissueExpressionData.mat')
load('.\TissueInitialModel.mat')
load('.\cons_ids.mat')
tissues = expressionData.Tissue;
Cmodel=removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],a)));
clear model

% fastcore
fn = '.\csm_models\StanDep\FASTCORE\';
load('.\SD_inputs\FASTCORE')
for i=1:54
    curr_core = core(a,i);
    model = fastcore(Cmodel,find(curr_core),1e-8);
    save([fn,tissues{i}],'model')
end
clear core

% gimme
fn = '.\csm_models\StanDep\GIMME\';
load('.\SD_inputs\GIMME')
for i=1:54
    curr_core = core(a,i);
    model = GIMME(Cmodel,curr_core,0);
    save([fn,tissues{i}],'model')
end
clear core

% imat
fn = '.\csm_models\StanDep\iMAT\';
load('.\SD_inputs\iMAT')
for i=1:54
    curr_core = core(a,i);
    model = iMAT(Cmodel,curr_core,0.2,0.7,1e-8,{});
    save([fn,tissues{i}],'model')
end
clear core


% mba
fn = '.\csm_models\StanDep\MBA\';
load('.\SD_inputs\MBA')
for i=1:54
    curr_H = H(a,i);
    curr_M = M(a,i);
    model = MBA(Cmodel,Cmodel.rxns(find(curr_M)),Cmodel.rxns(find(curr_H)),1e-8);
    save([fn,tissues{i}],'model')

end
clear H M

% init
fn = '.\csm_models\StanDep\INIT\';
load('.\SD_inputs\INIT')
for i=1:54
    curr_core = core(a,i);
    model = INIT(Cmodel,curr_core,1e-8);
    save([fn,tissues{i}],'model')
end
clear core

% mCADRE
fn = '.\csm_models\StanDep\mCADRE\';
load('.\SD_inputs\mCADRE')
for i=1:54
    curr_core = core(a,i);
    conf_scr = zeros(numel(a),1);
    model = mCADRE(Cmodel,curr_core,conf_scr,{},0,1/3,1e-8);
    save([fn,tissues{i}],'model')
end
clear core