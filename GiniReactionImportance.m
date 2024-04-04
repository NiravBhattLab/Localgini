function [RxnImp,Contexts] = GiniReactionImportance(geneExpression,model,MeM,ut,lt,ThS,coreRxn)

%%INPUT
%       geneExpression: matlab structure with fields
%                       .value : mRNA expression matrix with dimension N_genes*N_samples
%                       .genes : cell array with geneIDs in the same format as model.genes
%                       .context : cell array with names of the samples
%
%       model: COBRA model structure
%
%       MeM: Char variable indicating the model extraction method
%            has to be anyone of ['FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE']

%%OPTIONAL INPUTS        
%       ut: Scalar value [0 100] (upper threshold percentile) (default:100)
%         
%       lt: Scalar value [0 100] (lower threshold percentile) (default:0)
%         
%       ThS: Scalar variable that determines where thresholding has to be given
%              1 : Gene level ==>Thresholding to gene_expression followed by mapping to reactions
%              2 : Enzyme level ==>Thresholding to enzyme_expression followed by mapping to reactions
%              3 : Reaction level ==>Maps the gene_expression to reactions followed by thresholding
%              (default value: 1)
%         
%       coreRxn: Integer vector denoting indices of reactions for which 
%                higher importance has to be given manually (default:none)

%%OUTPUT
%       RxnImp: Reaction importance matrix with dimension N_genes*N_samples
%          FASTCORE: elements will be 1 or 0 defining core/non-core reactions
%          iMAT: elements will be 10 or 0 defining core/non-core reactions, all other reactions will carry 5
%          MBA: elements will be 2 or 1 defining High/Medium confidence reactions, all other reactions will carry 0
%          GIMME: elements are GIMME scores
%          INIT: elements are INIT weights
%          mCADRE: elements are ubiquity scores in the range of [0 1]

%%AUTHOR
%       Pavan Kumar S, BioSystems Engineering and Control (BiSECt) lab, IIT Madras


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



Contexts = geneExpression.context;

% getting expression values of genes available in the model
[idxa,~] = ismember(geneExpression.genes,model.genes);
model_exp_value = geneExpression.value(idxa,:);
model_exp_gene = geneExpression.genes(idxa);
parsedGPR = GPRparser(model);

if ThS==1 % thresholding at gene level
    
    temp_var = model_exp_value;
    temp_var(temp_var==0)=[];    
    LT = prctile(temp_var,lt); % lower threshold value
    UT = prctile(temp_var,ut); % upper threshold value   
    clear temp_var
    Loc_Gini = ginicoeff(model_exp_value); % gini coefficient based thresholding
    Loc_Gini(Loc_Gini>=UT)=UT;Loc_Gini(Loc_Gini<=LT)=LT;
    gene_exp = model_exp_value-repmat(Loc_Gini,1,numel(Contexts));
    gene_exp = gene_exp./std(gene_exp,[],2); % modified gene expression values    
    
    
    for i=1:numel(Contexts)
        rxn_exp(:,i) = selectGeneFromGPR(model, model_exp_gene, gene_exp(:,i), parsedGPR); % gene to reaction mapping
    end
    
elseif ThS==2 % thresholding at enzyme level
  
    for i=1:numel(Contexts)
        [enz_exp(:,i),enzymes,newparsedGPR] = geneToEnzymeExp(model_exp_value(:,i), model_exp_gene, parsedGPR); % gene to enzyme mapping
    end
    temp_var = enz_exp;
    temp_var(temp_var==0)=[];    
    LT = prctile(temp_var,lt); % lower threshold value
    UT = prctile(temp_var,ut); % upper threshold value   
    clear temp_var
    Loc_Gini = ginicoeff(enz_exp); % gini coefficient based thresholding
    Loc_Gini(Loc_Gini>=UT)=UT;Loc_Gini(Loc_Gini<=LT)=LT;
    enz_exp = enz_exp-repmat(Loc_Gini,1,numel(Contexts));
    enz_exp = enz_exp./std(enz_exp,[],2); % modified enzyme expression values    
       
    for i=1:numel(Contexts)
        rxn_exp(:,i) = enzymeToRxnExp(enz_exp(:,i),enzymes,newparsedGPR); % enzyme to reaction mapping
    end
    
elseif ThS==3 % thresholding at reaction level
    
    for i=1:numel(Contexts)
        rxn_exp(:,i) = selectGeneFromGPR(model, model_exp_gene, model_exp_value(:,i), parsedGPR); % gene to reaction mapping
    end
    temp_var = rxn_exp;
    temp_var(temp_var==0)=[];    
    LT = prctile(temp_var,lt); % lower threshold value
    UT = prctile(temp_var,ut); % upper threshold value   
    clear temp_var
    Loc_Gini = ginicoeff(rxn_exp); % gini coefficient based thresholding
    Loc_Gini(Loc_Gini>=UT)=UT;Loc_Gini(Loc_Gini<=LT)=LT;
    rxn_exp = rxn_exp-repmat(Loc_Gini,1,numel(Contexts));
    rxn_exp = rxn_exp./std(rxn_exp,[],2); % modified reaction expression values    
        
    
else
    error("error: ThS has to be one of [1,2,3]")
end

if ~isfield(model,'rxnGeneMat')
    model=buildRxnGeneMat(model);
end

if strcmp(MeM,'FASTCORE')
    
    RxnImp = rxn_exp;
    RxnImp(rxn_exp>=0)=1; RxnImp(rxn_exp<0)=0; RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=0;    
    RxnImp(coreRxn,:)=1;
    
elseif strcmp(MeM,'iMAT')
    
    RxnImp = rxn_exp;
    RxnImp(rxn_exp>=0)=10; RxnImp(rxn_exp<0)=0; RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=5;
    RxnImp(coreRxn,:)=10;
    
elseif strcmp(MeM,'MBA')
    
    RxnImp = rxn_exp;
    thr = prctile(RxnImp(RxnImp>=0),50,'all');
    RxnImp(rxn_exp>=thr)=2; RxnImp(rxn_exp<0)=0; RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=0;
    RxnImp(rxn_exp<thr & rxn_exp>=0)=1;
    RxnImp(coreRxn,:)=2;
    
elseif strcmp(MeM,'GIMME')
    
    eps=0.001;
    RxnImp = rxn_exp;
    RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=-eps;RxnImp(rxn_exp<-10)=-10;
    RxnImp(coreRxn,:)=max(RxnImp,[],'all');
    
elseif strcmp(MeM,'INIT')
    
    RxnImp = rxn_exp;
    RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=0;RxnImp(rxn_exp<-10)=-10;
    RxnImp(coreRxn,:)=max(RxnImp,[],'all');
    
elseif strcmp(MeM,'mCADRE')
    
    eps=0.001;
    RxnImp = rxn_exp;
    RxnImp(rxn_exp>=0)=1;RxnImp(find(sum(model.rxnGeneMat,2)==0),:)=1-eps;RxnImp(rxn_exp<-10)=-10;
    
    temp_var = RxnImp(RxnImp<0);    
    mxl=max(temp_var);mnl=min(temp_var);
    temp_var = (((temp_var)- mnl)/(mxl- mnl))*(1-eps);
    RxnImp(RxnImp<0) = temp_var;
    clear temp_var
    RxnImp(coreRxn,:)=1;
    
else
    error("error: MeM has to be one of ['FASTCORE','iMAT','MBA','GIMME','INIT','mCADRE']")
end
end

function [gcp]=ginicoeff(data)
    data_sort=sort(data,2);
    n=size(data,2);
    
    for i=1:size(data,1)
        x=data_sort(i,:);
        G_num=sum(((2*[1:n])-n-1).*x);
        G_den=sum(x)*n;
        gc(i)=(G_num/G_den)*100;
    end
    
    for i=1:numel(gc)
        gcp(i)=prctile(data(i,:),gc(i));
    end
    gcp=gcp';
end 

function [enz_exp,enzymes,newparsedGPR] = geneToEnzymeExp(geneExpression, geneNames, parsedGPR)
    enzymes=cell(0,0); % to store enzyme names
    enz_exp=[]; % to store the expression values of enzymes
    newparsedGPR=cell(size(parsedGPR));
    for i=1:numel(parsedGPR)
        current_rxn = parsedGPR{i};
        if ~isempty(current_rxn{1})
            newparsedGPR{i}=cell(size(current_rxn));
            for j=1:numel(current_rxn)
                currentArr = current_rxn{j};
                if length(currentArr)==1
                    newparsedGPR{i}{j}=currentArr;
                    if sum(ismember(enzymes,currentArr{1}))<1
                        ids = find(ismember(geneNames,currentArr{1}));
                        if ~isempty(ids)
                            currentExp = geneExpression(ids);
                            enzymes=[enzymes;currentArr{1}];
                            enz_exp=[enz_exp;currentExp];
                        end
                    end
                else
                    enz_name = strjoin(sort(currentArr),' and ');
                    newparsedGPR{i}{j} = cellstr(enz_name);
                    if sum(ismember(enzymes,enz_name))<1
                        ids=find(ismember(geneNames,currentArr));
                        if numel(ids)==numel(currentArr)
                            enzymes=[enzymes;enz_name];
                            currentExp = min(geneExpression(ids));
                            enz_exp=[enz_exp;currentExp];
                        end
                    end
                end
            end  
        else
            newparsedGPR{i}=current_rxn;
        end
        
    end   
end

function rxn_exp = enzymeToRxnExp(enz_exp,enzymes,parsedGPR)
    rxn_exp = -1*ones(length(parsedGPR),1);
    for i = 1:length(parsedGPR)
        current_rxn=parsedGPR{i};
        curArr=[];
        for j=1:length(current_rxn)
            if length(current_rxn{j})>=1
                ids = find(ismember(enzymes,current_rxn{j}));
                if ~isempty(ids) 
                     curArr= [curArr, enz_exp(ids)]; 
                end
            end
        end
        if ~isempty(curArr)
            rxn_exp(i)=max(curArr);
        end
    end
end
