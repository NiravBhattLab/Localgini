clear

%% for iHuman model
load('iHuman')
load('HK_rxns_iHuman')
hk_ids = find(ismember(model.rxns,HK_rxns));
load('core_rxns_cancer_data','LG_core_ihuman')
if sum(sum(LG_core_ihuman(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end
load('core_rxns_tissue_data','LG_core_ihuman')
if sum(sum(LG_core_ihuman(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end

%% for recon22 model
load('Recon2_2')
load('HK_rxns_Recon2_2')
hk_ids = find(ismember(model.rxns,HK_rxns));
load('core_rxns_cancer_data','LG_core_recon22')
if sum(sum(LG_core_recon22(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end
load('core_rxns_tissue_data','LG_core_recon22')
if sum(sum(LG_core_recon22(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end
%% for recon3d model
load('Recon3D')
load('HK_rxns_Recon3D')
hk_ids = find(ismember(model.rxns,HK_rxns));
load('core_rxns_cancer_data','LG_core_recon3d')
if sum(sum(LG_core_recon3d(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end
load('core_rxns_tissue_data','LG_core_recon3d')
if sum(sum(LG_core_recon3d(hk_ids,:),2)~=0)~=numel(hk_ids)
    fprintf('\n LG does not have all the HK rxns')
end