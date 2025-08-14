clear
cancer = load('core_rxns_cancer_data');
tissue = load('core_rxns_tissue_data');

c_ih = getOP_SScount(cancer.LG_core_ihuman,'iHuman');
c_r2 = getOP_SScount(cancer.LG_core_recon22,'Recon2_2');
c_r3d = getOP_SScount(cancer.LG_core_recon3d,'Recon3D');

t_ih = getOP_SScount(tissue.LG_core_ihuman,'iHuman');
t_r2 = getOP_SScount(tissue.LG_core_recon22,'Recon2_2');
t_r3d = getOP_SScount(tissue.LG_core_recon3d,'Recon3D');

function tbl = getOP_SScount(core,GEM)
load(['HK_rxns_',GEM])
load(GEM)
ids = ismember(model.rxns,HK_rxns);
core = core(ids,:);
HK_rxns = model.rxns(ids);
HK_rxnNames = model.rxnNames(ids);

SS = model.subSystems(ids);
if iscell(SS{1})
    for i=1:numel(SS)
        SS{i} = SS{i}{1};
    end
end
op_ids = ismember(SS,'Oxidative phosphorylation');
OP_count = sum(core(op_ids,:),2);
OP_rxnNames = HK_rxnNames(op_ids);
OP_rxns = HK_rxns(op_ids);
tbl = table(OP_count,OP_rxns,OP_rxnNames);
end
