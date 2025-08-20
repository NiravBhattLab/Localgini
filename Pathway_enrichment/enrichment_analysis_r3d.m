clear
core = load('../HK_reactions_analysis/core_rxns_tissue_data.mat');
lg_core = core.LG_core_recon3d; %%%%%%%
sd_core = core.SD_core_recon3d;%%%%%%%
lt2_core = core.LT2_core_recon3d;%%%%%%%
load('..\InputData\Recon3D.mat')%%%%%%%
load('../InputData/TissueExpressionData.mat') 
enrich_lg = getEnrichMat(lg_core,expressionData,model);
enrich_lt2 = getEnrichMat(lt2_core,expressionData,model);
enrich_sd = getEnrichMat(sd_core,expressionData,model);
save('enrichment_mats_recon3d','enrich_lg','enrich_lt2','enrich_sd')

function enrich_mat = getEnrichMat(core,expressionData,model)
    core = double(core); 
    path_tiss=readtable('./path_tiss_pairs.xlsx');
    %%%% modify based on the GEM %%%%%%%
    path_tiss.Pathway(strcmp(path_tiss.Pathway,'Arginine and Proline Metabolism'))={'Arginine and proline metabolism'};
    path_tiss.Pathway(strcmp(path_tiss.Pathway,'beta-Alanine metabolism'))={'Beta-Alanine metabolism'};
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



