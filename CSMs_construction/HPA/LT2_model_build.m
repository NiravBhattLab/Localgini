clear
% initCobraToolbox(false)
load('../../InputData/TissueExpressionData.mat')
load('.\TissueInitialModel.mat')
load('.\cons_ids.mat')
tissues = expressionData.Tissue;

Cmodel=removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),a)));
clear model

% fastcore
fn = 'csm_models\LocalT2\FASTCORE\';
load('LT2_inputs\FASTCORE')
for i=1:54
    curr_core = core(a,i);
    model = fastcore(Cmodel,find(curr_core),1e-8);
    save([fn,tissues{i}],'model')
end
clear core model

% gimme
fn = 'csm_models\LocalT2\GIMME\';
load('LT2_inputs\GIMME')
for i=1:54
    curr_core = core(a,i);
    model = GIMME(Cmodel,curr_core,5*log(2));
    save([fn,tissues{i}],'model')
end
clear core model

% imat
fn = 'csm_models\LocalT2\iMAT\';
load('LT2_inputs\GIMME')
for i=1:54
    curr_core = core(a,i);
    model = iMAT(Cmodel,curr_core,5*log(2),5*log(2),1e-8,{});
    save([fn,tissues{i}],'model')
end
clear core model

% mba
fn = 'csm_models\LocalT2\MBA\';
load('LT2_inputs\MBA')
for i=1:54
    curr_H = H(a,i);
    curr_M = M(a,i);
    model = MBA(Cmodel,Cmodel.rxns(curr_M),Cmodel.rxns(curr_H),1e-8);
    save([fn,tissues{i}],'model')
end
clear H M model

% init
fn = 'csm_models\LocalT2\INIT\';
load('LT2_inputs\INIT')
for i=1:54
    curr_core = core(a,i);
    model = INIT(Cmodel,curr_core,1e-8);
    save([fn,tissues{i}],'model')
end
clear core model

% mCADRE
fn = 'csm_models\LocalT2\mCADRE\';
load('LT2_inputs\GIMME')
for i=1:54
    curr_core = core(a,i);
    conf_scr = zeros(numel(a),1);
    model = mCADRE(Cmodel,curr_core,conf_scr,[],[],[],1e-8);
    save([fn,tissues{i}],'model')
end
clear core model
