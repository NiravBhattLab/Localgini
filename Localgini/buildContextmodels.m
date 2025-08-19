function [Models,RxnImp] = buildContextmodels(geneExpression,model,MeM,contexts,ut,lt,ThS,coreRxn,filename,cons_mod_rxn_id,varargin)

%%INPUT
%       geneExpression: matlab structure with fields
%                       .value : mRNA gene expression matrix with dimension N_genes*N_samples
%                       .genes : cell array with geneIDs in the same format as model.genes
%                       .context : cell array with names of the samples
%
%       model: COBRA model structure
%
%       MeM: Char variable indicating the model extraction method
%            has to be anyone of ['FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE']
%
%       contexts: contexts for which models has to be built 
%               (has to be in the same format as geneExpression.context)

%%OPTIONAL INPUTS        
%       ut: Scalar value [0 100] (upper threshold percentile) (default:100)
%         
%       lt: Scalar value [0 100] (lower threshold percentile) (default:0)
%         
%       ThS: Scalar variable that determines when thresholding has to be given
%              1 : Gene level ==>Thresholding to gene_expression followed by mapping to reactions
%              2 : Enzyme level ==>Maps to enzyme_expression followed by thresholding to enzyme_expression followed by mapping to reactions
%              3 : Reaction level ==>Maps the gene_expression to reactions followed by thresholding to reaction_expression
%              (default value: 1)
%         
%       coreRxn: Integer vector denoting indices of reactions for which 
%                higher importance has to be given manually (default:none)
%
%       filename: Folder name to save the generated models in .mat format
%
%       cons_mod_rxn_id: Integer vector indicating the reaction ids of consistent reactions
%
%       varargin: Optional inputs required for different MeMs
%                 If GIMME ==> tolerance and objective_fraction
%                 If mCADRE ==> tolerance and confidence scores
%                 For all other MeMs ==> tolerance

%%OUTPUT
%       Models: Context specific model in COBRA structure
%       RxnImp: Reaction important scores derived for distinct MeMs by Localgini method

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% Setting the default values
if ~exist('ut','var') || isempty(ut)
    ut =100;
end

if ~exist('lt','var') || isempty(lt)
    lt =0;
end

if ~exist('ThS','var') || isempty(ThS)
    ThS =1;
end

if ~exist('coreRxn','var') || isempty(coreRxn)
    coreRxn =[];
end

if nargin>10
    tol = varargin{1};
else
    tol = 1e-4;
end

% flux consistent model
if isempty(cons_mod_rxn_id)
    cons_mod_rxn_id = [1:numel(model.rxns)];
else
    cons_mod_rxn_id = sort(unique(cons_mod_rxn_id));
end

cons_mod = removeRxns(model,model.rxns(setdiff([1:numel(model.rxns)],cons_mod_rxn_id)));

if strcmp(MeM,'FASTCORE')
    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);    
    
    for i=1:numel(contexts)
       current_core = find(RxnImp(:,i));   
       model=fastcore(cons_mod,current_core,tol);
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end     
    
elseif strcmp(MeM,'iMAT')
    
    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);    
    
    for i=1:numel(contexts)
       current_core = RxnImp(:,i);
       
       model=iMAT(cons_mod,current_core,1,9,tol,{});
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end
     
elseif strcmp(MeM,'MBA')
    
    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);    
    
    for i=1:numel(contexts)
       current_core = RxnImp(:,i);   
       model=MBA(cons_mod,cons_mod.rxns(current_core==1),cons_mod.rxns(current_core==2),tol);
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end
    
elseif strcmp(MeM,'GIMME')
    
    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);    
    
    if nargin>11
        obj_frac = varargin{2};
    else
        obj_frac = 0.9;
    end
    
    for i=1:numel(contexts)
       current_core = RxnImp(:,i);   
       model=GIMME(cons_mod,current_core,0,obj_frac);
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end    
    
elseif strcmp(MeM,'INIT')

    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);    

    for i=1:numel(contexts)
       current_core = RxnImp(:,i);   
       model = INIT(cons_mod,current_core,tol);
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end
    
elseif strcmp(MeM,'mCADRE')
    
    [RxnImp,~] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn);
    col = ismember(geneExpression.context,contexts);
    % picking the consistent reactions for the context specified
    RxnImp = RxnImp(cons_mod_rxn_id,col);  
    if nargin>11
        conf_scr = varargin{2}; % confidence scores
    else
        conf_scr = zeros(size(RxnImp,1),1); % confidence scores
    end
    
    
    for i=1:numel(contexts)
       current_core = RxnImp(:,i);   
       model = mCADRE(cons_mod,current_core,conf_scr,[],[],[],tol);
       if exist('filename','var')
            % saving the model
            fn = strcat(filename,contexts{i});
            save(fn,'model');
       end
       Models{i} = model;
    end
    
end
