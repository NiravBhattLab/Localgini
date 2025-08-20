clear
load('../InputData/CancerExpressionData.mat') %%%%%%
load('../CSMs_construction/CancerCellLines/generic_models/786O.mat') %%%%%
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

thrs = {'LocalT2','Localgini','StanDep'};
thrs2 = {'LT2','LocalGiniModels','StanDep'};
mems = {'FASTCORE','GIMME','iMAT','INIT','MBA','mCADRE'};

for t =1:3
    thr = thrs{t};
    sims=[];
    for c=1:44
        cl = cell_lines{c};
        sim_mat = zeros(6,numel(model.rxns));
        for m=1:6
            mem=mems{m};
            load(['../CSMs_construction/CancerCellLines/csm_models/',thr,'/',...
                mem,'/',cl,'.mat'])
            sim_mat(m,:) = ismember(model.rxns,m1.rxns);
        end
        C = nchoosek(1:6,2);
        for i=1:15
            sim = Jacc_sim(sim_mat(C(i,1),:),sim_mat(C(i,2),:));            
            sims=[sims,sim];
        end
    end
    save(['./cancer_models/JC_',thrs2{t}],'sims')
end
            