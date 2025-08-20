clear
load('../InputData/TissueExpressionData.mat')
load('../HK_reactions_analysis/core_rxns_tissue_data.mat')
m1 = load('..\CSMs_construction\HPA\TissueInitialModel.mat');
m1=m1.model;
load('..\CSMs_construction\HPA\cons_ids.mat')
m1=removeRxns(m1,m1.rxns(setdiff(1:numel(m1.rxns),a)));

tissues = expressionData.Tissue;

thrs = {'LocalT2','StanDep','Localgini'};
thrs2 = {'LT2','SD','LG'};

mems= {'FASTCORE','GIMME','iMAT','INIT','MBA','mCADRE'};

for t =1:3
    thr = thrs{t};
    eval(['core=',thrs2{t},'_core_recon22;']);
    for m=1:6
        mem=mems{m};
        frac_added_mat=[];
        for i=1:54
            core_rxns = m1.rxns(find(core(:,i)));
            load(['..\CSMs_construction\HPA\csm_models\',thr,'\',mem,'\',tissues{i},'.mat'])
            frac_added = numel(setdiff(model.rxns,core_rxns))/numel(model.rxns);
            frac_added_mat=[frac_added_mat frac_added];
            clear model
        end
        eval(['FA_',thrs2{t},'_',mems{m},'=frac_added_mat;'])        
    end
end
FA_LGgen_FC = FA_LG_FASTCORE;
FA_LGgen_GIMME = FA_LG_GIMME;
FA_LGgen_iMAT = FA_LG_iMAT;
FA_LGgen_INIT = FA_LG_INIT;
FA_LGgen_MBA = FA_LG_MBA;
FA_LGgen_mCADRE = FA_LG_mCADRE;

FA_LT2_FC = FA_LT2_FASTCORE;
FA_SD_FC = FA_SD_FASTCORE;


pval_fc_lt2 = ranksum(FA_LGgen_FC,FA_LT2_FC,'tail','left');
pval_gimme_lt2 = ranksum(FA_LGgen_GIMME,FA_LT2_GIMME,'tail','left');
pval_iMAT_lt2 = ranksum(FA_LGgen_iMAT,FA_LT2_iMAT,'tail','left');
pval_MBA_lt2 = ranksum(FA_LGgen_MBA,FA_LT2_MBA,'tail','left');
pval_INIT_lt2 = ranksum(FA_LGgen_INIT,FA_LT2_INIT,'tail','left');
pval_mCADRE_lt2 = ranksum(FA_LGgen_mCADRE,FA_LT2_mCADRE,'tail','left');

pval_fc_sd = ranksum(FA_LGgen_FC,FA_SD_FC,'tail','left');
pval_gimme_sd = ranksum(FA_LGgen_GIMME,FA_SD_GIMME,'tail','left');
pval_iMAT_sd = ranksum(FA_LGgen_iMAT,FA_SD_iMAT,'tail','left');
pval_MBA_sd = ranksum(FA_LGgen_MBA,FA_SD_MBA,'tail','left');
pval_INIT_sd = ranksum(FA_LGgen_INIT,FA_SD_INIT,'tail','left');
pval_mCADRE_sd = ranksum(FA_LGgen_mCADRE,FA_SD_mCADRE,'tail','left');

save('HPA/pvals_self_cons','pval_fc_lt2','pval_gimme_lt2','pval_iMAT_lt2','pval_MBA_lt2',...
    'pval_INIT_lt2','pval_mCADRE_lt2','pval_fc_sd','pval_gimme_sd','pval_iMAT_sd',...
    'pval_MBA_sd','pval_INIT_sd','pval_mCADRE_sd')
save('HPA/self_cons_plots','FA_LGgen_FC','FA_LGgen_GIMME','FA_LGgen_iMAT','FA_LGgen_INIT','FA_LGgen_MBA','FA_LGgen_mCADRE',...
    'FA_LT2_FC','FA_LT2_GIMME','FA_LT2_iMAT','FA_LT2_INIT','FA_LT2_MBA','FA_LT2_mCADRE',...
    'FA_SD_FC','FA_SD_GIMME','FA_SD_iMAT','FA_SD_INIT','FA_SD_MBA','FA_SD_mCADRE')