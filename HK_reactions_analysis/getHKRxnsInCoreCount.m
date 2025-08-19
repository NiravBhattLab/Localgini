clear
% for cancer data
load('core_rxns_cancer_data.mat')
load('../InputData/CancerExpressionData.mat')
contexts = expressionData.Tissue;
thresholds = {'LG','LT2','SD'};
GEMs = {'ihuman','recon22','recon3d'};

for t=1:numel(thresholds)
    t = thresholds{t};
    for g=1:numel(GEMs)
        g = GEMs{g};
        if strcmp(g, 'ihuman')
            temp = 'iHuman';
        elseif strcmp(g, 'recon22')
            temp = 'Recon2_2';
        elseif strcmp(g, 'recon3d')
            temp = 'Recon3D';
        end
        eval(['core =',t,'_core_',g,';']) 
        eval(['load("HK_rxns_' temp '")'])
        eval(['load("../InputData/',temp,'")'])
        ids = ismember(model.rxns,HK_rxns);
        hk_count = sum(core(ids,:),1);
        eval(['nHK_',t,'_',g,'=hk_count;'])
    end
end
save('HK_in_core_cancer','contexts','nHK_LG_ihuman','nHK_LG_recon22','nHK_LG_recon3d','nHK_LT2_ihuman','nHK_LT2_recon22',...
'nHK_LT2_recon3d','nHK_SD_ihuman','nHK_SD_recon22','nHK_SD_recon3d')

% for tissue data
load('core_rxns_tissue_data.mat')
load('../InputData/TissueExpressionData.mat')
contexts = expressionData.Tissue;
thresholds = {'LG','LT2','SD'};
GEMs = {'ihuman','recon22','recon3d'};

for t=1:numel(thresholds)
    t = thresholds{t};
    for g=1:numel(GEMs)
        g = GEMs{g};
        if strcmp(g, 'ihuman')
            temp = 'iHuman';
        elseif strcmp(g, 'recon22')
            temp = 'Recon2_2';
        elseif strcmp(g, 'recon3d')
            temp = 'Recon3D';
        end
        eval(['core =',t,'_core_',g,';']) 
        eval(['load("HK_rxns_' temp '")'])
        eval(['load("../InputData/',temp,'")'])
        ids = ismember(model.rxns,HK_rxns);
        hk_count = sum(core(ids,:),1);
        eval(['nHK_',t,'_',g,'=hk_count;'])
    end
end
save('HK_in_core_tissue','contexts','nHK_LG_ihuman','nHK_LG_recon22','nHK_LG_recon3d','nHK_LT2_ihuman','nHK_LT2_recon22',...
'nHK_LT2_recon3d','nHK_SD_ihuman','nHK_SD_recon22','nHK_SD_recon3d')