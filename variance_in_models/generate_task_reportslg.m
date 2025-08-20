%% this code generates metabolic task report for all the built models
clear
initCobraToolbox(false)
load('../InputData/CancerExpressionData.mat') 
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';


mems = {'FASTCORE','GIMME','mCADRE','MBA','iMAT','INIT'};

for m = 1:numel(mems)
    mem = mems{m};
    f = ['../CSMs_construction/CancerCellLines/csm_models/Localgini/',mem,'/task_rep/'];
    for c = 1:44
        load(['../CSMs_construction/CancerCellLines/csm_models/Localgini/',mem,'/',cell_lines{c}]);
        if sum(ismember(m1.rxns,{'sink_citr(c)','DM_citr_L[c]'}))>1
            m1 = removeRxns(m1,{'sink_citr(c)'});
        end
        if sum(ismember(m1.rxns,{'DM_anth'  ,'DM_anth[c]'}))>1
            m1 = removeRxns(m1,{'DM_anth'});
        end
        m1.lb(ismember(m1.rxns,{'biomass_reaction','DM_atp_c_'}))=0;
        task_rep = checkMetabolicTasks(m1,'TASKS.xlsx');
        save([f,'/',cell_lines{c}],'task_rep')
    end
end
            