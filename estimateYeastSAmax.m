load('GEM-yeast.mat');
model_split = splitRevRxns(model);
model = model_split;

model = setYeastMedia(model);
model = blockYeastRxns(model);
model = changeRxnBounds(model, 'r_4046', 0, 'l');
model = changeRxnBounds(model, 'r_4046', 1000, 'u');
model = changeObjective(model, 'r_4046');

% import molecular weight
[num,txt,~] = xlsread('UniProt_Yeast.xlsx');
geneUniProt = txt(2:end,1);
MWUniProt = num;
clear num txt;

%% 
[num,txt,~] = xlsread('Yeast_fluxes_proteomics.xlsx','fluxes');
exRxnList = txt(2:end,1);
exFluxes = num;
expList = txt(1,3:end);
clear num txt;

SA_estimated = zeros(length(model.rxns),length(expList));
flux_simulated = zeros(length(model.rxns),length(expList));
for i = 1:length(expList) % whether exclude DiBartolomeo data
    display([num2str(i),'/',num2str(length(expList))]);
    
    expID = expList{i};
    sheetName = expID(1:strfind(expID,'_')-1);
    [num,txt,~] = xlsread('Yeast_fluxes_proteomics.xlsx',sheetName);
    head = txt(1,2:end);
    proteinList = txt(2:end,1);
    abundList = num(:,ismember(head,expID));
    
    exFluxes_tmp = exFluxes(:,i)';
    exRxnList_tmp = exRxnList';
    
    idxnan = isnan(exFluxes_tmp);
    exFluxes_tmp = exFluxes_tmp(~idxnan);
    exRxnList_tmp = exRxnList_tmp(~idxnan);
    
    model_tmp = model;
    
    % change protein content for N-limited chemostats
    if contains(expID,'Yu_CN30')
        model_tmp = scaleBioMass(model_tmp,'protein',0.3665,'carbohydrate');
    elseif contains(expID,'Yu_CN50') || contains(expID,'Yu_CN115')
        model_tmp = scaleBioMass(model_tmp,'protein',0.2635,'carbohydrate');
    end
    
    % set max NGAM for batch conditions
    if contains(expID,'DiBartolomeo')
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0, 'l');
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0.7, 'u');
    end
    
    [fluxes, deviation] = searchFeasibleSol(model_tmp,exRxnList_tmp,exFluxes_tmp,'max',0.01);
    % remove large deviation
    if deviation > 0.2
        fluxes(:,1) = 0;
    end
    
    fluxes(abs(fluxes)<1e-3) = 0; %﻿the absolute flux value should surpass 0.001
    
    flux_simulated(:,i) = fluxes;
    
    for j = 1:length(fluxes)
        if fluxes(j) ~= 0 && ~ismember(model_tmp.grRules(j),'')
            gr_tmp = model_tmp.grRules{j};
            gr_tmp = strrep(gr_tmp,'(','');
            gr_tmp = strrep(gr_tmp,')','');
            gr_tmp = strrep(gr_tmp,' and ',' ');
            gr_tmp = strrep(gr_tmp,' or ',' ');
            gr_tmp = strtrim(gr_tmp);
            gr_tmp = strsplit(gr_tmp);
            gr_tmp = unique(gr_tmp)';
            
            [a,~] = ismember(gr_tmp,proteinList);
            if any(a)
                gr_tmp = gr_tmp(a);
                if strcmp(sheetName,'Lahtvee')
                    [~,abund_idx] = ismember(gr_tmp,proteinList);
                    abund_tmp = abundList(abund_idx);
                    [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                    mw_tmp = MWUniProt(mw_idx);
                    mol_gCDW = abund_tmp*1e12/6.02e23;
                    mass = sum(mw_tmp.*mol_gCDW)*1000; %mg/gCDW
                    sa = fluxes(j)*1000/60/mass; %umol/mg/min
                elseif strcmp(sheetName,'Yu')
                    [~,abund_idx] = ismember(gr_tmp,proteinList);
                    abund_tmp = abundList(abund_idx);
                    [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                    mw_tmp = MWUniProt(mw_idx);
                    mol_gCDW = abund_tmp*1e3/1e15;
                    mass = sum(mw_tmp.*mol_gCDW)*1000; %mg/gCDW
                    sa = fluxes(j)*1000/60/mass; %umol/mg/min
                elseif strcmp(sheetName,'DiBartolomeo')
                    [~,abund_idx] = ismember(gr_tmp,proteinList);
                    abund_tmp = abundList(abund_idx);
                    mass = sum(abund_tmp)*1000; %mg/gCDW
                    sa = fluxes(j)*1000/60/mass; %umol/mg/min
                end
                SA_estimated(j,i) = sa;
            end
        end
    end
end
SA_estimated(SA_estimated == inf) = 0;
SA_max = max(SA_estimated,[],2);

YeastSAmax.rxn = model.rxns(SA_max~=0);
YeastSAmax.samax = SA_max(SA_max~=0);
YeastSAmax.grRules = model.grRules(SA_max~=0);
YeastSAmax.condition = cell(length(YeastSAmax.rxn),1);

for i = 1:length(YeastSAmax.rxn)
    all_sa = SA_estimated(ismember(model.rxns,YeastSAmax.rxn(i)),:);
    YeastSAmax.condition(i,1) = expList(all_sa == max(all_sa));
end


save('YeastSAmax.mat','YeastSAmax');
