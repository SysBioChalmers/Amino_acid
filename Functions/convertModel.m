%% convertModel 
function model = convertModel(model,enzymedata)
% integrate protein cost into a split model
% protein costs will be added into enzymatic reactions as metabolites

cutoff = 0;

model = addMetabolite(model,'protein_cost[c]');
met_idx = ismember(model.mets,'protein_cost[c]');

for i = 1:length(enzymedata.rxn)
    rxnid = enzymedata.rxn(i);
    rxn_idx = ismember(model.rxns,rxnid);
    if enzymedata.kcat_conf(i) >= cutoff
        coeff = -enzymedata.minMW(i)/enzymedata.kcat(i);
        model.S(met_idx,rxn_idx) = coeff;
    end
end

model = addReaction(model,'prot_cost_exchange',...
                    'metaboliteList',{'protein_cost[c]'},...
                    'stoichCoeffList',1000,...
                    'reversible',false,...
                    'lowerBound',0,...
                    'upperBound',1000);






