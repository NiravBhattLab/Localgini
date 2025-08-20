function [JS,frac_added] = SelfConsistency(Con_gen_model,Context_model,core_rxns)


%%INPUT
%       Con_gen_model: Consistent genome scale model (COBRA model stucture)
%
%       Context_model: Context specific model (COBRA model structure)
%
%       core_rxns: Consistent core reactions in same format as
%                  Context_model.rxns

%%OUTPUT
%       JS: Jaccard similarity between core reaction list and reactions from context-specific model
%                 
%       frac_added: Fractional contribution of reactions by MeMs

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% Jaccard similarity between core reactions and model reactions
core_rxn_vec =ismember(Con_gen_model.rxns,core_rxns);
mod_rxn_vec =ismember(Con_gen_model.rxns,Context_model.rxns);
JS = Jacc_sim(core_rxn_vec,mod_rxn_vec);

% Fraction of reactions added by MeM in the context specific model
frac_added = numel(setdiff(Context_model.rxns,core_rxns))/numel(Context_model.rxns);
end

function similarity = Jacc_sim(A,B)
 similarity = sum(A.*B)/(sum(A+B) - sum(A.*B)); % Jaccard similarity of A,B
end
