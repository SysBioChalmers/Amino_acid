%% updateYeastEnzyme
function enzymedata = updateYeastEnzyme(enzymedata,model)

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

% update kcat values based on manual check
[num, txt, ~] = xlsread('manualYeast.xlsx','kcat');
rxnlist = txt(2:end,1);
kcatmanual = num(:,1);
kcatconf = num(:,2);
for i = 1:length(rxnlist)
    if ismember(rxnlist(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,rxnlist(i));
        if enzymedata.kcat_conf(idx_tmp) < kcatconf(i)
            enzymedata.kcat(idx_tmp) = 3600*kcatmanual(i);
            enzymedata.kcat_conf(idx_tmp) = kcatconf(i);
        elseif enzymedata.kcat_conf(idx_tmp) == kcatconf(i) && enzymedata.kcat(idx_tmp) < 3600*kcatmanual(i)
            enzymedata.kcat(idx_tmp) = 3600*kcatmanual(i);
        end
    end
end

% update kcat values from https://github.com/SysBioChalmers/Yeast_kapp
% from collected in vitro kcat
load('../Yeast_kapp/kcat.mat');
for i = 1:length(kcat.rxn)
    if ismember(kcat.rxn(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,kcat.rxn(i));
        if enzymedata.kcat_conf(idx_tmp) < 4
            enzymedata.kcat(idx_tmp) = 3600*kcat.value(i);
            enzymedata.kcat_conf(idx_tmp) = 4;
        elseif enzymedata.kcat_conf(idx_tmp) == 4 && enzymedata.kcat(idx_tmp) < 3600*kcat.value(i)
            enzymedata.kcat(idx_tmp) = 3600*kcat.value(i);
        end
    end
end

% add protein cost
enzymedata.minPC = enzymedata.minMW./enzymedata.kcat;

load('../Yeast_kapp/kapp.mat');
for i = 1:length(kapp.rxn)
    if ismember(kapp.rxn(i),enzymedata.rxn)
        idx_tmp = ismember(enzymedata.rxn,kapp.rxn(i));
        if enzymedata.kcat_conf(idx_tmp) < 5 && max(kapp.values(i,:)) > 0
            enzymedata.kcat(idx_tmp) = max(kapp.values(i,:))*3600;
            enzymedata.minPC(idx_tmp) = enzymedata.minMW(idx_tmp)/(max(kapp.values(i,:))*3600);
            enzymedata.kcat_conf(idx_tmp) = 5;
        elseif enzymedata.kcat_conf(idx_tmp) == 5 && enzymedata.kcat(idx_tmp) < max(kapp.values(i,:))*3600
            enzymedata.kcat(idx_tmp) = max(kapp.values(i,:))*3600;
            enzymedata.minPC(idx_tmp) = enzymedata.minMW(idx_tmp)/(max(kapp.values(i,:))*3600);
        end
    end
end



