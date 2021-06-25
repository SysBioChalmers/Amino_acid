function enzymedata = updateCondSpecYeastEnzyme(enzymedata,condID)

load('../Yeast_kapp/kapp.mat');
kappraw = kapp.values(:,ismember(kapp.condition,condID));
idx = kappraw ~= 0;
rxnlist = kapp.rxn(idx);
kapplist = kappraw(idx);
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        if enzymedata.kcat_conf(idx_tmp) < 6
            enzymedata.kcat(idx_tmp) = kapplist(i)*3600;
            enzymedata.kcat_conf(idx_tmp) = 6;
            enzymedata.minPC(idx_tmp) = enzymedata.minMW(idx_tmp)/(kapplist(i)*3600);
        elseif enzymedata.kcat_conf(idx_tmp) == 6 && enzymedata.kcat(idx_tmp) < kapplist(i)*3600
            enzymedata.kcat(idx_tmp) = kapplist(i)*3600;
            enzymedata.minPC(idx_tmp) = enzymedata.minMW(idx_tmp)/(kapplist(i)*3600);
        end
    end
end

% use mu-dependent kmax
kmaxothers = max(kapp.values(:, kapp.condGR <= kapp.condGR(ismember(kapp.condition,condID))),[],2);
idxnew = (kmaxothers ~= 0) & ~idx;
rxnlist = kapp.rxn(idxnew);
kapplist = kmaxothers(idxnew);
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        enzymedata.kcat(idx_tmp) = kapplist(i)*3600;
        enzymedata.kcat_conf(idx_tmp) = 5.5;
        enzymedata.minPC(idx_tmp) = enzymedata.minMW(idx_tmp)/(kapplist(i)*3600);
    end
end







