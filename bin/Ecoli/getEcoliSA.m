load('GEM-ecoli-split.mat');

[num, txt, ~] = xlsread('sa_ecoli_Davidi.xlsx');
rxnlist = txt(2:end,1);
values = num;
values(isnan(values)) = 0;
EcoliSA.condition = txt(1,2:end);
EcoliSA.values = zeros(0,length(EcoliSA.condition));
EcoliSA.protein = cell(0,1);
EcoliSA.rxn = cell(0,1);

for i = 1:length(rxnlist)
    rxnid_tmp = rxnlist(i);
    rxnid_tmp = strrep(rxnid_tmp,'_reverse','_rvs');
    if ismember(strcat(rxnid_tmp,'_fwd'),model_split.rxns)
        rxnid_tmp = strcat(rxnid_tmp,'_fwd');
        EcoliSA.rxn = [EcoliSA.rxn;rxnid_tmp];
        EcoliSA.values = [EcoliSA.values;values(i,:)];
        EcoliSA.protein = [EcoliSA.protein;model_split.grRules(ismember(model_split.rxns,rxnid_tmp))];
    elseif ismember(rxnid_tmp,model_split.rxns)
        EcoliSA.rxn = [EcoliSA.rxn;rxnid_tmp];
        EcoliSA.values = [EcoliSA.values;values(i,:)];
        EcoliSA.protein = [EcoliSA.protein;model_split.grRules(ismember(model_split.rxns,rxnid_tmp))];
    end
end

save('EcoliSA.mat','EcoliSA');





