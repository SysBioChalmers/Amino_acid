% Timing: ~ 2400 s
tic;

load('../Yeast_kapp/kappEstimation/GEM-yeast-split.mat');
rxnlist = model_split.rxns;
grlist = model_split.grRules;
rowlist = cell(1,0);
valuelist = zeros(length(rxnlist),0);

% import molecular weight
[num,txt,~] = xlsread('../Yeast_kapp/kappEstimation/Data/UniProt.xlsx');
geneUniProt = txt(2:end,1);
MWUniProt = num;
clear num txt;

cd ../Yeast_kapp/kappEstimation/Fluxes/;
file = dir('*.mat');
for i = 1:length(file)
    display([num2str(i) '/' num2str(length(file))]);
    filename = file(i).name;
    load(filename);
    if Fluxes.deviation < 0.1
        condid = strrep(filename,'Fluxes_','');
        condid = strrep(condid,'.mat','');
        rowlist = [rowlist {condid}];
        
        sheetName = condid(1:strfind(condid,'_')-1);
        [num,txt,~] = xlsread('../Data/ProteomicsFlux.xlsx',sheetName);
        % remove proteins with extremely low abundance (< 1 percentile)
        numtmp = num(:);
        numtmp = numtmp(numtmp ~= 0);
        num(num < quantile(numtmp,0.01)) = 0;
        
        head = txt(1,2:end);
        proteinList = txt(2:end,1);
        abundList = num(:,ismember(head,condid));
        
        idxnonzero = abundList ~= 0;
        abundList = abundList(idxnonzero);
        proteinList = proteinList(idxnonzero);
        
        fluxes = Fluxes.pFBA;
        fluxes(abs(fluxes)<1e-5) = 0; %﻿the absolute flux value should surpass 0.00001
        salist_tmp = zeros(length(rxnlist),1);
        for j = 1:length(rxnlist)
            flux_tmp = fluxes(ismember(Fluxes.model.rxns,rxnlist(j)));
            
            if flux_tmp > 0 && ~ismember(grlist(j),'')
                tot_flux_tmp = Fluxes.pFBA(ismember(Fluxes.model.rxns,'r_tot_flux'));
                model_tmp = Fluxes.model;
                model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*0.9999,'l');
                model_tmp = changeRxnBounds(model_tmp,'r_tot_flux',tot_flux_tmp*1.0001,'u');
                model_tmp = changeObjective(model_tmp,rxnlist{j});
                solmin = optimizeCbModel(model_tmp,'min');
                
                % remove data when min of FVA is zero
                if solmin.f > 0
                    gr_tmp = grlist{j};
                    gr_tmp = strrep(gr_tmp,'(','');
                    gr_tmp = strrep(gr_tmp,')','');
                    gr_tmp = strrep(gr_tmp,' or ',' ');
                    gr_tmp = strrep(gr_tmp,' and ',' ');
                    gr_tmp = strtrim(gr_tmp);
                    gr_tmp = strsplit(gr_tmp);
                    gr_tmp = unique(gr_tmp)';

                    [a,~] = ismember(gr_tmp,proteinList);
                    if (any(a) && ~contains(grlist{j},' and ')) || (all(a) && contains(grlist{j},' and '))
                        gr_tmp = gr_tmp(a);
                        if strcmp(sheetName,'Lahtvee2017')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                            mw_tmp = MWUniProt(mw_idx);
                            mol_gCDW = abund_tmp*1e12/6.02e23;
                            mass = sum(mw_tmp.*mol_gCDW)*1000; %mg/gCDW
                            sa_tmp = flux_tmp*1000/60/mass; %umol/mg/min
                        elseif strcmp(sheetName,'Yu2020') || strcmp(sheetName,'Yu2021')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                            mw_tmp = MWUniProt(mw_idx);
                            mol_gCDW = abund_tmp*1e3/1e15;
                            mass = sum(mw_tmp.*mol_gCDW)*1000; %mg/gCDW
                            sa_tmp = flux_tmp*1000/60/mass; %umol/mg/min
                        elseif strcmp(sheetName,'DiBartolomeo2020')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mass = sum(abund_tmp)*1000; %mg/gCDW
                            sa_tmp = flux_tmp*1000/60/mass; %umol/mg/min
                        end
                        salist_tmp(j,1) = sa_tmp;
                        
                    elseif ismember(rxnlist(j),{'r_0226' 'r_0438' 'r_0439' 'r_0961' 'r_1021_fwd'})
                        [num,txt,~] = xlsread('YeastComplex.xlsx',rxnlist{j});
                        list = txt(2:end,1);
                        s = num(:,1);
                        mw = num(:,2);
                        mol = zeros(length(list),1);
                        
                        gr_tmp = intersect(gr_tmp(a),list);
                        if strcmp(sheetName,'Lahtvee2017')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mol_gCDW = abund_tmp*1e12/6.02e23;
                        elseif strcmp(sheetName,'Yu2020') || strcmp(sheetName,'Yu2021')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            mol_gCDW = abund_tmp*1e3/1e15;
                        elseif strcmp(sheetName,'DiBartolomeo2020')
                            [~,abund_idx] = ismember(gr_tmp,proteinList);
                            abund_tmp = abundList(abund_idx);
                            [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                            mw_tmp = MWUniProt(mw_idx);
                            mol_gCDW = abund_tmp./mw_tmp;
                        end
                        
                        [~,s_idx] = ismember(gr_tmp,list);
                        s_tmp = s(s_idx);
                        mol_undetected = mean(mol_gCDW./s_tmp); % select mean/min/median?
                        for l = 1:length(mol)
                            if ismember(list(l),gr_tmp)
                                mol(l) = mol_gCDW(ismember(gr_tmp,list(l)));
                            else
                                mol(l) = mol_undetected * s(l);
                            end
                        end
                        
                        mass = sum(mw.*mol)*1000; %mg/gCDW
                        sa_tmp = flux_tmp*1000/60/mass; %umol/mg/min
                        salist_tmp(j,1) = sa_tmp;
                    end
                end
            end
        end
        valuelist = [valuelist salist_tmp];
    end
end
valuelist(valuelist == inf) = 0;

conditions = strrep(rowlist,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
conditions = unique(conditions);

YeastSA = struct();
YeastSA.values = zeros(length(rxnlist),length(conditions));
YeastSA.rxn = rxnlist;
YeastSA.protein = grlist;
YeastSA.condition = conditions;

for i = 1:length(conditions)
    idx = contains(rowlist,conditions(i));
    values = valuelist(:,idx);
    values(values == 0) = nan;
    if size(values,2) > 1
        for j = 1:size(values,1)
            cov = std(values(j,:),'omitnan')/mean(values(j,:),'omitnan');
            if cov > 0.5 || isnan(cov) % remove data with very high coefficient of variation
                YeastSA.values(j,i) = nan;
            else
                YeastSA.values(j,i) = max(values(j,:));
            end
        end
    elseif size(values,2) == 1
        YeastSA.values(:,i) = values;
    end
end
YeastSA.values(isnan(YeastSA.values)) = 0;

idx1 = any(YeastSA.values,2); % remain rxns that have sa at >= 1 condition
YeastSA.rxn = YeastSA.rxn(idx1);
YeastSA.protein = YeastSA.protein(idx1);
YeastSA.values = YeastSA.values(idx1,:);

cd ../../../Amino_acid/; 
save('YeastSA.mat','YeastSA');

toc;

