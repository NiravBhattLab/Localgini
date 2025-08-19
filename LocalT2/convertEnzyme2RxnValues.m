function [rxnData,enzSelected] = convertEnzyme2RxnValues(enzymeData,model)

n = size(enzymeData.value,2);
specialist = find(~cellfun(@iscell,enzymeData.rxns));
generalist = find(cellfun(@iscell,enzymeData.rxns));

k = 0;
generalistActivity = [];
generalistRxns = [];
g = [];
for i=1:length(generalist)
    t = length(enzymeData.rxns{generalist(i)});
    g = [g;repmat(generalist(i),t,1)];
    generalistActivity(end+1:end+t,:) = repmat(enzymeData.value(generalist(i),:),t,1);
    generalistRxns = [generalistRxns;enzymeData.rxns{generalist(i)}];
end
s = specialist;
specialistActivity = enzymeData.value(specialist,:);
specialistRxns = enzymeData.rxns(specialist);
rxnActivity = [generalistActivity;specialistActivity];
e = [g;s];
enzs = enzymeData.enzyme(e);
rxns = [generalistRxns;specialistRxns];
rxns = cellstr(rxns);
rxnList = unique(rxns);
modelActivity = repmat(-2,length(model.rxns),n);
for i=1:length(rxnList)
    [modelActivity(strcmp(model.rxns,rxnList{i}),:),ie] = max(rxnActivity(strcmp(rxns,rxnList{i}),:),[],1);
    renz = enzs(strcmp(rxns,rxnList{i}));
    enzSelected(strcmp(model.rxns,rxnList{i}),:) = renz(ie);
end
rxnData.value = modelActivity;
rxnData.rxns = model.rxns;
rxnData.Tissue = enzymeData.Tissue;