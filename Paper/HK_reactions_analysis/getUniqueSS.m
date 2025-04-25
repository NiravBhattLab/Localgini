clear
%% For cancer data
load('core_rxns_cancer_data')
% for the case of ihuman
[c_ss_ih,c_rxns_ih] = getSusystemsLocalGini(LG_core_ihuman,LT2_core_ihuman,SD_core_ihuman,'iHuman');
% for the case of Recon2.2
[c_ss_r2,c_rxns_r2] = getSusystemsLocalGini(LG_core_recon22,LT2_core_recon22,SD_core_recon22,'Recon2_2');
% for the case of Recon3D
[c_ss_r3,c_rxns_r3] = getSusystemsLocalGini(LG_core_recon3d,LT2_core_recon3d,SD_core_recon3d,'Recon3D');

%% For tissue data
load('core_rxns_tissue_data')
% for the case of ihuman
[t_ss_ih,t_rxns_ih] = getSusystemsLocalGini(LG_core_ihuman,LT2_core_ihuman,SD_core_ihuman,'iHuman');
% for the case of Recon2.2
[t_ss_r2,t_rxns_r2] = getSusystemsLocalGini(LG_core_recon22,LT2_core_recon22,SD_core_recon22,'Recon2_2');
% for the case of Recon3D
[t_ss_r3,t_rxns_r3] = getSusystemsLocalGini(LG_core_recon3d,LT2_core_recon3d,SD_core_recon3d,'Recon3D');


function [ss,rxns] = getSusystemsLocalGini(lg,lt2,sd,modelName)
load(modelName)
load(['HK_rxns_',modelName])
hk_ids = find(ismember(model.rxns,HK_rxns));
lg_ = find(sum(lg,2));
lt2_ = find(sum(lt2,2));
sd_ = find(sum(sd,2));

lg_unique = setdiff(lg_,sd_);
lg_unique = setdiff(lg_unique,lt2_);
ss = model.subSystems(lg_unique);
rxns = model.rxns(lg_unique);
end

