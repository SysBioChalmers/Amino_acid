function enzymedata = updateCondSpecEcoliEnzyme(enzymedata,condID)
load('EcoliSA.mat');
% idx = all(EcoliSA.values,2);
% EcoliSA.protein = EcoliSA.protein(idx);
% EcoliSA.rxn = EcoliSA.rxn(idx);
% EcoliSA.values = EcoliSA.values(idx,:);

saraw = EcoliSA.values(:,ismember(EcoliSA.condition,condID));

idx = saraw ~= 0;

rxnlist = EcoliSA.rxn(idx);
salist = saraw(idx);

for i = 1:length(rxnlist)
%     if ismember(rxnlist(i),enzymedata.rxn) && ~contains(EcoliSA.protein(i),' and ')
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        if enzymedata.kcat_conf(idx_tmp) <= 6
            enzymedata.minPC(idx_tmp) = 1000/60/salist(i);
            enzymedata.kcat_conf(idx_tmp) = 6;
        elseif enzymedata.kcat_conf(idx_tmp) == 6 && enzymedata.minPC(idx_tmp) > 1000/60/salist(i)
            enzymedata.minPC(idx_tmp) = 1000/60/salist(i);
        end
    end
end