function [accuracy, pvals, sgd]  = GeneEssentialityAcc(model,crispr_scores,cutoff)
%%INPUT
%       model: context-specific model COBRA model structure
%
%       crispr_scores: matlab structure with fields
%                       .value : CRISPR-CAS9 loss of function screens of the corresponding context
%                       .genes : cell array with geneIDs in the same format as model.genes
%       cutoff: growth rate ratio between mutant strain and wild type below which the genes are considered to be essential

%%OUTPUT
%       accuracy: (Number of genes predicted to be essential with -ve CRISPR scores)/(Total number of predicted essential genes)
%       pval: pvalue returned after wilcoxon rank-sum test with H1 as essential genes have less scores than non-essential genes
%       sgd: genes predicted to be essential
%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

gene_idx = ismember(crispr_scores.genes,model.genes);
crispr_scores.genes = crispr_scores.genes(gene_idx);
crispr_scores.value = crispr_scores.value(gene_idx);
sgd = fastSL_sg(model,cutoff);
ess_gene_scores = crispr_scores.value(ismember(crispr_scores.genes,sgd));
non_ess_gene_scores = crispr_scores.value(~ismember(crispr_scores.genes,sgd));
accuracy = sum(ess_gene_scores<0)/numel(ess_gene_scores);
pvals = ranksum(ess_gene_scores,non_ess_gene_scores,'tail','left');

end
