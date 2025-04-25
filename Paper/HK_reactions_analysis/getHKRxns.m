function HK_rxns = getHKRxns(GEM)
% To get the house keeping (HK) reactions from the given input genome-scale
% metabolic model (GEM)
%
% USAGE:
%   HK_rxns = getHKRxns(GEM)
%
% INPUTS:
%   GEM:     Name of the model ('Recon2_2','Recon3D','iHuman')
%
% OUTPUTS:
%   HK_rxns: House keeping reactions in the model


%% Dataset loading
HK_gene = importdata('./hkg_sym.txt');% Housekeeping genes data
load(GEM); % GEM

if strcmp(GEM,'Recon2_2')
    field = 'genes';
    filename = './gecko_id_to_entrez.txt';% gene symbol to entrz id
    opts = detectImportOptions(filename);
    symtoentrez = readtable(filename,opts);
    HK_gene = symtoentrez.Var2(ismember(symtoentrez.Var1,HK_gene));
    HK_gene = arrayfun(@num2str,HK_gene,'UniformOutput',false);
elseif strcmp(GEM,'Recon3D')
    field = 'geneisrefseq_nameID';
    model = buildRxnGeneMat(model);
elseif strcmp(GEM,'iHuman')
    field = 'geneShortNames';
end

met_HK_genes = intersect(model.(field),HK_gene);
ids = ismember(model.(field),met_HK_genes);
HK_rxns_ids = sum(model.rxnGeneMat(:,ids),2)>=1;
HK_rxns = model.rxns(HK_rxns_ids);
save(['HK_rxns_',GEM],'HK_rxns')
end
