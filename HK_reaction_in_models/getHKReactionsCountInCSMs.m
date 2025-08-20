clear
% for cancer cell line data
load('..\InputData\CancerExpressionData.mat')
load('..\HK_reactions_analysis\HK_rxns_Recon2_2.mat')
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';
thrs = {'LocalT2','StanDep','Localgini'};
thrs2 = {'LT2','StanDep','LocalGiniModels'};
mems = {'FASTCORE','iMAT','INIT','GIMME','mCADRE','MBA'};

for t1 =1:numel(thrs)
    t = thrs{t1};
    for m = 1:numel(mems)
        m = mems{m};
        fracs = [];
        for i=1:numel(cell_lines)
            load(['../CSMS_construction/CancerCellLines/csm_models/',t,'/',m,'/',cell_lines{i},'.mat'])
            frac = sum(ismember(m1.rxns,HK_rxns))/numel(HK_rxns);
            fracs(i) = frac;
        end
        eval([thrs2{t1},'_',m,'=fracs;'])
    end
end

save('CancerCellLines/HK_model.mat','LT2_FASTCORE','LT2_iMAT','LT2_INIT','LT2_GIMME','LT2_mCADRE','LT2_MBA',...
'StanDep_FASTCORE','StanDep_iMAT','StanDep_INIT','StanDep_GIMME','StanDep_mCADRE','StanDep_MBA',...
'LocalGiniModels_FASTCORE','LocalGiniModels_iMAT','LocalGiniModels_INIT','LocalGiniModels_GIMME','LocalGiniModels_mCADRE','LocalGiniModels_MBA')

% for Tissue data
load('..\InputData\TissueExpressionData.mat')
load('..\HK_reactions_analysis\HK_rxns_Recon2_2.mat')
tissue_types = expressionData.Tissue;

thrs = {'LocalT2','StanDep','Localgini'};
thrs2 = {'LT2','StanDep','LocalGiniModels'};
mems = {'FASTCORE','iMAT','INIT','GIMME','mCADRE','MBA'};

for t1 =1:numel(thrs)
    t = thrs{t1};
    for m = 1:numel(mems)
        m = mems{m};
        fracs = [];
        for i=1:numel(tissue_types)
            load(['../CSMS_construction/HPA/csm_models/',t,'/',m,'/',tissue_types{i},'.mat'])
            frac = sum(ismember(model.rxns,HK_rxns))/numel(HK_rxns);
            fracs(i) = frac;
        end
        eval([thrs2{t1},'_',m,'=fracs;'])
    end
end

save('HPA/HK_model.mat','LT2_FASTCORE','LT2_iMAT','LT2_INIT','LT2_GIMME','LT2_mCADRE','LT2_MBA',...
'StanDep_FASTCORE','StanDep_iMAT','StanDep_INIT','StanDep_GIMME','StanDep_mCADRE','StanDep_MBA',...
'LocalGiniModels_FASTCORE','LocalGiniModels_iMAT','LocalGiniModels_INIT','LocalGiniModels_GIMME','LocalGiniModels_mCADRE','LocalGiniModels_MBA')