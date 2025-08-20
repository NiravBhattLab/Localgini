clear
load('../InputData/TissueExpressionData.mat') %%%%%%
load('../CSMs_construction/HPA/TissueInitialModel.mat') %%%%%
gem = model;
cell_lines = expressionData.Tissue;

thrs = {'LocalT2','Localgini','StanDep'};
thrs2 = {'LT2','LocalGiniModels','StanDep'};
mems = {'FASTCORE','GIMME','iMAT','INIT','MBA','mCADRE'};

for t =1:3
    thr = thrs{t};
    sims=[];
    for c=1:54
        cl = cell_lines{c};
        sim_mat = zeros(6,numel(gem.rxns));
        for m=1:6
            mem=mems{m};
            load(['../CSMs_construction/HPA/csm_models/',thr,'/',...
                mem,'/',cl,'.mat'])
            sim_mat(m,:) = ismember(gem.rxns,model.rxns);
        end
        C = nchoosek(1:6,2);
        for i=1:15
            sim = Jacc_sim(sim_mat(C(i,1),:),sim_mat(C(i,2),:));            
            sims=[sims,sim];
        end
    end
    save(['./tissue_models/JC_',thrs2{t}],'sims')
end
            