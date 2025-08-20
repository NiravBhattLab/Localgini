clear
core = load('../HK_reactions_analysis/core_rxns_tissue_data.mat');
lg_core = load('..\CSMs_construction\HPA\LG_inputs\FASTCORE.mat');
lg_core = lg_core.core_rxn;
sd_core = load('..\CSMs_construction\HPA\SD_inputs\FASTCORE.mat');
sd_core = sd_core.core;
lt2_core = load('..\CSMs_construction\HPA\LT2_inputs\FASTCORE.mat');
lt2_core = lt2_core.core;

load('../CSMs_construction/HPA/TissueInitialModel.mat')
load('../InputData/TissueExpressionData.mat') 
enrich_lg = getEnrichMat(lg_core,expressionData,model);
enrich_lt2 = getEnrichMat(lt2_core,expressionData,model);
enrich_sd = getEnrichMat(sd_core,expressionData,model);
save('enrichment_mats_recon22','enrich_lg','enrich_lt2','enrich_sd')

function enrich_mat = getEnrichMat(core,expressionData,model)
    core = double(core); 
    path_tiss=readtable('./path_tiss_pairs.xlsx');
    pathways=unique(path_tiss.Pathway);
    tissues=intersect(path_tiss.Tissue,expressionData.Tissue);
    n_tissue = numel(tissues);
    n_pathway= numel(pathways);
    enrich_mat=zeros(n_tissue,n_pathway);
    M=numel(model.rxns); % total reactions in the model
    for tiss=1:n_tissue
        id = find(ismember(expressionData.Tissue,tissues{tiss}));
        N=sum(core(:,id)); % number of core reactions 
        for path=1:n_pathway
            K=sum(strcmp(model.subSystems,pathways{path})); % number of reactions in the pathway
            x=sum(ismember(find(core(:,id)),find(strcmp(model.subSystems,pathways{path}))));
            p=hygecdf(x,M,K,N,'upper');
            enrich_mat(tiss,path)=p<0.05;
        end
    end
end



