clear
%% for cancer data
data = load('core_rxns_cancer_data');
cancer_ih = getHKcount(data.LG_core_ihuman,'iHuman');
cancer_r2 = getHKcount(data.LG_core_recon22,'Recon2_2');
cancer_r3d = getHKcount(data.LG_core_recon3d,'Recon3D');

writetable(cancer_ih, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Cancer_ih');
writetable(cancer_r2, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Cancer_r2');
writetable(cancer_r3d, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Cancer_r3d');

%% for tissue data
data = load('core_rxns_tissue_data');
tissue_ih = getHKcount(data.LG_core_ihuman,'iHuman');
tissue_r2 = getHKcount(data.LG_core_recon22,'Recon2_2');
tissue_r3d = getHKcount(data.LG_core_recon3d,'Recon3D');

writetable(tissue_ih, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Tissue_ih');
writetable(tissue_r2, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Tissue_r2');
writetable(tissue_r3d, 'HK_tissue_Countdata.xlsx', 'Sheet', 'Tissue_r3d');

function tbl = getHKcount(core,GEM)
load(['HK_rxns_',GEM])
load(GEM)
ids = ismember(model.rxns,HK_rxns);
HKcount = sum(core(ids,:),2);
SS = model.subSystems(ids);
rxns = model.rxns(ids);
tbl = table(HKcount,SS,rxns);
end

% function plotBySSGrouping(count,GEM)
% load(['HK_rxns_',GEM])
% load(GEM)
% ids = ismember(model.rxns,HK_rxns);
% HK_rxns = model.rxns(ids);
% SS = model.subSystems(ids);
% uSS = unique(SS);
% uSS = uSS(~cellfun(@isempty,uSS));
% i=1;c_id=[true,false];
% figure()
% cs={'b','r'};
% for s =1:numel(uSS)
%     ss_ids= find(ismember(SS,uSS{s}));
%     c = cs{c_id};
%     c_id =~c_id;
%     for r =1:numel(ss_ids)
%         hold on;
%         bar(i,count(ss_ids(r)),c,'LineWidth',0.000001)
%         i=1+i;
%     end
% end
% end