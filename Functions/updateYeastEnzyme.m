%% updateYeastEnzyme
function enzymedata = updateYeastEnzyme(enzymedata,model,expID)

[num, txt, ~] = xlsread('UniProt_Yeast.xlsx');
unip_gene = txt(2:end,1);
unip_mw = num;

enzymedata.minMW = zeros(length(enzymedata.rxn),1);

% add minimal enzyme molecular weight
for i = 1:length(enzymedata.rxn)
    genes_tmp = model.grRules{ismember(model.rxns,enzymedata.rxn(i))};
    if ~isempty(genes_tmp)
        iso_tmp = strsplit(genes_tmp,' or ');
        iso_tmp = strrep(iso_tmp,'(','');
        iso_tmp = strrep(iso_tmp,')','');
        iso_tmp = strtrim(iso_tmp);
        mw_tmp = zeros(1,length(iso_tmp));
        for j = 1:length(iso_tmp)
            iso_tmp_tmp = iso_tmp{j};
            if ~contains(iso_tmp_tmp,'and')
                if ismember(iso_tmp_tmp,unip_gene)
                    mw_tmp(1,j) = unip_mw(ismember(unip_gene,iso_tmp_tmp));
                end
            else
                subunit_tmp = strsplit(iso_tmp_tmp,' and ');
                subunit_tmp = strrep(subunit_tmp,'(','');
                subunit_tmp = strrep(subunit_tmp,')','');
                subunit_tmp = strtrim(subunit_tmp);
                [~,b] = ismember(subunit_tmp,unip_gene);
                mw_tmp(1,j) = sum(unip_mw(b));
            end
        end
        if any(mw_tmp)
            enzymedata.minMW(i,1) = min(mw_tmp(mw_tmp ~= 0));
        end
    end
end

% update molecular weight for complexes
[num, txt, ~] = xlsread('manualYeast.xlsx','mw');
rxnlist = txt(2:end,1);
mwmanual = num;
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        enzymedata.minMW(ismember(enzymedata.rxn,rxnlist(i))) = mwmanual(i);
    end
end


% update kcat values based on estimated max specific activity 
% only greater one adopted compared to the data
% load('YeastSA.mat');
% idx = ~contains(YeastSA.grRules,'and');
% YeastSA.grRules = YeastSA.grRules(idx);
% YeastSA.rxn = YeastSA.rxn(idx);
% YeastSA.samax = YeastSA.samax(idx);
% 
% for i = 1:length(YeastSA.rxn)
%     rxnid_tmp = YeastSA.rxn(i);
%     minMW = enzymedata.minMW(i);
%     samax = YeastSA.samax(i);
%     kmax_tmp = samax*minMW*60/1000;
%     
%     idx_tmp = ismember(enzymedata.rxn,rxnid_tmp);
%     if enzymedata.kcat_conf(idx_tmp) < 4
%         enzymedata.kcat(idx_tmp) = kmax_tmp;
%         enzymedata.kcat_conf(idx_tmp) = 5;
%     else
%         if kmax_tmp > enzymedata.kcat(idx_tmp)
%             enzymedata.kcat(idx_tmp) = kmax_tmp;
%             enzymedata.kcat_conf(idx_tmp) = 5;
%         end
%     end
% end

load(['YeastSA_' expID '.mat']);
idx = ~contains(YeastSA.grRules,'and');
YeastSA.grRules = YeastSA.grRules(idx);
YeastSA.rxn = YeastSA.rxn(idx);
YeastSA.sa = YeastSA.sa(idx);
for i = 1:length(YeastSA.rxn)
    rxnid_tmp = YeastSA.rxn(i);
    minMW = enzymedata.minMW(i);
    samax = YeastSA.sa(i);
    kmax_tmp = samax*minMW*60/1000;
    
    idx_tmp = ismember(enzymedata.rxn,rxnid_tmp);
    enzymedata.kcat(idx_tmp) = kmax_tmp;
    enzymedata.kcat_conf(idx_tmp) = 5;
end


% update kcat values based on manual check
[num, txt, ~] = xlsread('manualYeast.xlsx','kcat');
rxnlist = txt(2:end,1);
kcatmanual = num(:,1);
kcatconf = num(:,2);
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        if enzymedata.kcat_conf(idx_tmp) <= kcatconf(i)
            enzymedata.kcat(idx_tmp) = 3600*kcatmanual(i);
            enzymedata.kcat_conf(idx_tmp) = kcatconf(i);
        end
    end
end











