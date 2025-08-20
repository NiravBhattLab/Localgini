clear
load('../InputData/CancerExpressionData.mat')
[~,~,raw] = xlsread('hall_mark.xlsx');
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

thrs = {'LocalT2','StanDep','Localgini'};
mems = {'FASTCORE','GIMME','iMAT','INIT','MBA','mCADRE'};

hall_mark_gen = raw(3:199,2);
hall_mark_gen = cellfun(@num2str,hall_mark_gen,'UniformOutput',false);


for t =1:3
    thr = thrs{t};
    for m = 1:6
        mem =mems{m} ;
        eval([thr,'_',mem,'=[];'])
        for i=1:44  
            load(['../CSMs_construction/CancerCellLines/csm_models/',thr,'/',mem,'/',cell_lines{i},'.mat'])
            m1 = removeUnusedGenes(m1);
            acc_per = sum(ismember(hall_mark_gen,m1.genes))/numel(hall_mark_gen)*100;
            eval([thr,'_',mem,'=[',thr,'_',mem,' acc_per','];'])
        end
    end
end
save('per_hall_mark_rec','LocalGiniModels_FASTCORE','LocalGiniModels_GIMME','LocalGiniModels_iMAT','LocalGiniModels_INIT',...
    'LocalGiniModels_MBA','LocalGiniModels_mCADRE','LT2_FASTCORE','LT2_GIMME','LT2_iMAT','LT2_INIT',...
    'LT2_MBA','LT2_mCADRE','StanDep_FASTCORE','StanDep_GIMME','StanDep_iMAT','StanDep_INIT',...
    'StanDep_MBA','StanDep_mCADRE')
