clear
% initCobraToolbox(false)
load('../../InputData/CancerExpressionData.mat')
cell_lines = expressionData.Tissue;
cell_lines = regexprep(cell_lines,'_','','emptymatch');
cell_lines{strcmp(cell_lines,'x786O')}='786O';
cell_lines{strcmp(cell_lines,'NIHOVCAR3')}='OVCAR3';

thr = 'StanDep'; %%%% only change here
mems = {'FASTCORE','GIMME','iMAT','INIT','MBA','mCADRE'};

pca_mat =zeros(44*6,210);
count =1;
for m=1:6
    mem = mems{m};
    for i=1:44
        load(['../../CSMs_construction/CancerCellLines/csm_models/',thr,'/',...
                mem,'/task_rep/',cell_lines{i},'.mat'])
        pca_mat(count,:)=ismember(task_rep(1:end-2,5),'true');
        count =count+1;
        clear task_rep
    end
end

% removing columns with all reactions or no reactions present
pca_mat = pca_mat(:,sum(pca_mat,1)~=0);
pca_mat = pca_mat(:,sum(pca_mat,1)~=44*6);
% centering the pca_mat to 0 
pca_mat = pca_mat - repmat(mean(pca_mat,2),1,size(pca_mat,2));

[coeff,score,latent] = pca(pca_mat);


% calculating explained variance in first and second principal component for mem
vars_pc1 = []; vars_pc2 = []; vars_pc3 = []; % explained variance in all possibe orderings of the mem categories
pc1 = score(:,1);pc2 = score(:,2);pc3 = score(:,3);
all_poss_order = perms(1:6);
for i =1:size(all_poss_order,1)
    curr_order = all_poss_order(i,:);
    var_1 = [ones(44,1)*curr_order(1);ones(44,1)*curr_order(2);ones(44,1)*curr_order(3);ones(44,1)*curr_order(4);ones(44,1)*curr_order(5);ones(44,1)*curr_order(6)];
    co1=corrcoef(var_1,pc1);co2=corrcoef(var_1,pc2);co3=corrcoef(var_1,pc3);
    vars_pc1 = [vars_pc1;co1(2)];
    vars_pc2 = [vars_pc2;co2(2)];
    vars_pc3 = [vars_pc3;co3(2)];
end

explained_var_in_pc1_mem = max(vars_pc1.*vars_pc1*100);
explained_var_in_pc2_mem = max(vars_pc2.*vars_pc2*100);
explained_var_in_pc3_mem = max(vars_pc3.*vars_pc3*100);



% calculating explained variance in first and second principal component for cell lines
vars_pc1 = []; vars_pc2 = []; vars_pc3 = [];% explained variance in random 100000 orderings of the celllines
pc1 = score(:,1);pc2 = score(:,2);pc3 = score(:,3);
for i =1:10000
    curr_order = randperm(44);
    var_1 = repmat(curr_order,1,6);
    co1=corrcoef(var_1,pc1);co2=corrcoef(var_1,pc2);co3=corrcoef(var_1,pc3);
    vars_pc1 = [vars_pc1;co1(2)];
    vars_pc2 = [vars_pc2;co2(2)];
    vars_pc3 = [vars_pc3;co3(2)];
end

explained_var_in_pc1_cl = max(vars_pc1.*vars_pc1*100);
explained_var_in_pc2_cl = max(vars_pc2.*vars_pc2*100);
explained_var_in_pc3_cl = max(vars_pc3.*vars_pc3*100);



% calculating explained variance in first and second principal component for cancer type
[~,~,raw] = xlsread('cellline_cancer.xlsx');
cancer_type = raw(2:end,2); cell_line_map = raw(2:end,1);
for i = 1:44
    cancer_type_ordered{i,1} = cancer_type{find(ismember(cell_line_map,upper(cell_lines{i})))};
end
uni_cancer_type =unique(cancer_type_ordered);
vars_pc1 = []; vars_pc2 = [];vars_pc3 = []; % explained variance in random 100000 orderings of the cancer type
pc1 = score(:,1);pc2 = score(:,2);pc3 = score(:,3);
for i =1:10000
    curr_order_0 = randperm(numel(uni_cancer_type));
    curr_order = zeros(44,1);
    for j = 1:numel(uni_cancer_type)
       curr_order(ismember(cancer_type_ordered,uni_cancer_type{j})) =curr_order_0(j); 
    end
    var_1 = repmat(curr_order,1,6);
    co1=corrcoef(var_1,pc1);co2=corrcoef(var_1,pc2);co3=corrcoef(var_1,pc3);
    vars_pc1 = [vars_pc1;co1(2)];
    vars_pc2 = [vars_pc2;co2(2)];
    vars_pc3 = [vars_pc3;co3(2)];
end

explained_var_in_pc1_ct = max(vars_pc1.*vars_pc1*100);
explained_var_in_pc2_ct = max(vars_pc2.*vars_pc2*100);
explained_var_in_pc3_ct = max(vars_pc3.*vars_pc3*100);



save('all_mem/SD_var_met_task','explained_var_in_pc1_mem','explained_var_in_pc2_mem','explained_var_in_pc3_mem','explained_var_in_pc1_cl','explained_var_in_pc2_cl','explained_var_in_pc3_cl','explained_var_in_pc1_ct','explained_var_in_pc2_ct','explained_var_in_pc3_ct')