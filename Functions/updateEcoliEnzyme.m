%% updateEcoliEnzyme
function enzymedata = updateEcoliEnzyme(enzymedata,model)

[num, txt, ~] = xlsread('UniProt_Ecoli.xlsx');
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
[num, txt, ~] = xlsread('manualEcoli.xlsx','mw');
rxnlist = txt(2:end,1);
mwmanual = num;
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        enzymedata.minMW(ismember(enzymedata.rxn,rxnlist(i))) = mwmanual(i);
    end
end

% update kcat values based on max specific activity (DOI: 10.1073/pnas.1514240113)
% only greater one adopted compared to the data
[num, txt, ~] = xlsread('manualEcoli.xlsx','sa');
rxnlist = txt(2:end,1);
samanual = num;
for i = 1:length(rxnlist)
    rxnid_tmp = rxnlist(i);
    rxnid_tmp = strrep(rxnid_tmp,'_reverse','_rvs');
    sa_tmp = samanual(i);
    if ismember(rxnid_tmp,enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnid_tmp);
        mw_tmp = enzymedata.minMW(idx_tmp);
        kcandidate = 60*sa_tmp*mw_tmp/1000; % unit: /h
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = kcandidate;
            enzymedata.kcat_conf(idx_tmp) = 5;
        else
            if kcandidate > enzymedata.kcat(idx_tmp)
                enzymedata.kcat(idx_tmp) = kcandidate;
                enzymedata.kcat_conf(idx_tmp) = 5;
            end
        end
    elseif ismember(strcat(rxnid_tmp,'_fwd'),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,strcat(rxnid_tmp,'_fwd'));
        mw_tmp = enzymedata.minMW(idx_tmp);
        kcandidate = 60*sa_tmp*mw_tmp/1000; % unit: /h
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = kcandidate;
            enzymedata.kcat_conf(idx_tmp) = 5;
        else
            if kcandidate > enzymedata.kcat(idx_tmp)
                enzymedata.kcat(idx_tmp) = kcandidate;
                enzymedata.kcat_conf(idx_tmp) = 5;
            end
        end
    end
end

% update kcat values based on kmax (DOI: 10.1101/767996), only greater one
% adopted compared to the data with conf > 3
[num, txt, ~] = xlsread('manualEcoli.xlsx','kmax');
rxnlist = txt(2:end,1);
kmaxmanual = num;
for i = 1:length(rxnlist)
    rxnid_tmp = rxnlist(i);
    rxnid_tmp = strrep(rxnid_tmp,'_b','_rvs');
    kmax_tmp = kmaxmanual(i)*3600; % unit: /h
    if ismember(rxnid_tmp,enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnid_tmp);
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = kmax_tmp;
            enzymedata.kcat_conf(idx_tmp) = 5;
        else
            if kmax_tmp > enzymedata.kcat(idx_tmp)
                enzymedata.kcat(idx_tmp) = kmax_tmp;
                enzymedata.kcat_conf(idx_tmp) = 5;
            end
        end
    elseif ismember(strcat(rxnid_tmp,'_fwd'),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,strcat(rxnid_tmp,'_fwd'));
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = kmax_tmp;
            enzymedata.kcat_conf(idx_tmp) = 5;
        else
            if kmax_tmp > enzymedata.kcat(idx_tmp)
                enzymedata.kcat(idx_tmp) = kmax_tmp;
                enzymedata.kcat_conf(idx_tmp) = 5;
            end
        end
    end
end

% update kcat values based on manual check
[num, txt, ~] = xlsread('manualEcoli.xlsx','kcat');
rxnlist = txt(2:end,1);
kcatmanual = num(:,1);
kcatconf = num(:,2);
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = 3600*kcatmanual(i);
            enzymedata.kcat_conf(idx_tmp) = kcatconf(i);
        else
            if 3600*kcatmanual(i) > enzymedata.kcat(idx_tmp)
                enzymedata.kcat(idx_tmp) = 3600*kcatmanual(i);
                enzymedata.kcat_conf(idx_tmp) = kcatconf(i);
            end
        end
    end
end











