%% calculateProteinCost 
function [tot,pclist] = calculateProteinCost(enzymedata,rxnlist,sol,kcat_cutoff)

pclist = zeros(length(rxnlist),1);
for i = 1:length(rxnlist)
    kcat_conf = enzymedata.kcat_conf(ismember(enzymedata.rxn,rxnlist(i)));
    if kcat_conf >= kcat_cutoff
        kcat = enzymedata.kcat(ismember(enzymedata.rxn,rxnlist(i)));
        mw = enzymedata.minMW(ismember(enzymedata.rxn,rxnlist(i)));
        flux = sol(i);
        if flux >= 0
            pclist(i,1) = flux*mw/kcat/1000;
        end
    end
end
tot = sum(pclist);
