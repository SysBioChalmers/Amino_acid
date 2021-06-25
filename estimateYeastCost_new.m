% for each condition just calculate the condition-specific flux changes and
% protein costs.

% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% cd Models/;
% save('GEM-yeast-split.mat','model_split');
% cd ../;

load('GEM-yeast-split.mat');

enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);

step = 0.0001;
mu = 0.4;
substrateExID = 'r_1714'; % glucose

% kmax for estimating cost

model_max = convertModel(model_split,enzymedata);
model_max = setYeastMedia(model_max);
cost_kmax = calculateYeastCost(model_max,substrateExID,mu,step);
cd YeastCost/;
save('cost_kmax.mat','cost_kmax');
cd ../;


% condition-specific kapp for estimating cost

load('../Yeast_kapp/kapp.mat');
cond = kapp.condition(~contains(kapp.condition,'etoh','IgnoreCase',true));

for k = 1:length(cond)
    display([num2str(k),'/',num2str(length(cond))]);
    enzymedata_cond = enzymedata;
    enzymedata_cond = updateCondSpecYeastEnzyme(enzymedata_cond,cond(k));
    
    model = convertModel(model_split,enzymedata_cond);
    
    model = setYeastMedia(model);
    
    cost_kapp = calculateYeastCost(model,substrateExID,mu,step);
    
    
    cd YeastCost/;
    save(['cost_kapp_' cond{k} '.mat'],'cost_kapp');
    cd ../;
end



