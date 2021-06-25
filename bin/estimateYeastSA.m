% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% save('GEM-yeast-split.mat','model_split');

load('GEM-yeast-split.mat');
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

for i = 1:length(expList)
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
    elseif contains(expID,'std_010')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.053444 0.050594 0.045962]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.35967375 0.377942893 0.336456714]),'carbohydrate');
    elseif contains(expID,'N30_010')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.022625 0.02601 0.027257]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.209333929 0.221513357 0.245872214]),'carbohydrate');
    elseif contains(expID,'N30_013')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.026366 0.025831 0.025297]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.3037245 0.256909821 0.266425]),'carbohydrate');
    elseif contains(expID,'N30_018')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.039014 0.038658 0.037945]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.289261429 0.282029893 0.288119607]),'carbohydrate');
    elseif contains(expID,'N30_030')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.047031 0.045606 0.052019]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.422473929 0.421332107 0.3996375]),'carbohydrate');
    elseif contains(expID,'N30_035')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.073219 0.070546 0.088717]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.484154674 0.489921983 0.508864255]),'carbohydrate');
    elseif contains(expID,'Gln_glc')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.071259 0.058789 0.064133]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.497834143 0.583851357 0.59412775]),'carbohydrate');
    elseif contains(expID,'Gln_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.026366 0.022268 0.017815]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.205527857 0.352822821 0.361957393]),'carbohydrate');
    elseif contains(expID,'Phe_std')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.019596 0.020309 0.018171]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.425138179 0.520289964 0.452161286]),'carbohydrate');
    elseif contains(expID,'Phe_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.023337 0.024941 0.024406]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.417906643 0.378704107 0.398876286]),'carbohydrate');
    elseif contains(expID,'Ile_std')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.029216 0.034739 0.030641]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.320851821 0.418667857 0.331508821]),'carbohydrate');
    elseif contains(expID,'Ile_N30')
        model_tmp = scaleBioMass(model_tmp,'RNA',mean([0.024762 0.021912 0.024762]),'carbohydrate');
        model_tmp = scaleBioMass(model_tmp,'protein',mean([0.439981857 0.448355214 0.42361575]),'carbohydrate');
    end
    
    % set max NGAM for batch conditions
    if contains(expID,'DiBartolomeo')
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0, 'l');
        model_tmp = changeRxnBounds(model_tmp, 'r_4046', 0.7, 'u');
    end
    
    [fluxes, deviation] = searchFeasibleSol(model_tmp,exRxnList_tmp,exFluxes_tmp,'max',0.001);
    
    display(['deviation: ',num2str(deviation)]);
    
    fluxes(abs(fluxes)<1e-5) = 0; %﻿the absolute flux value should surpass 0.00001
    
    SA_estimated = zeros(length(model.rxns),1);
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
                elseif strcmp(sheetName,'Yu2021')
                    [~,abund_idx] = ismember(gr_tmp,proteinList);
                    abund_tmp = abundList(abund_idx);
                    [~,mw_idx] = ismember(gr_tmp,geneUniProt);
                    mw_tmp = MWUniProt(mw_idx);
                    mol_gCDW = abund_tmp*1e3/1e15;
                    mass = sum(mw_tmp.*mol_gCDW)*1000; %mg/gCDW
                    sa = fluxes(j)*1000/60/mass; %umol/mg/min
                end
                SA_estimated(j,1) = sa;
            end
        end
    end
    SA_estimated(SA_estimated == inf) = 0;
    YeastSA.deviation = deviation;
    YeastSA.rxn = model.rxns(SA_estimated~=0);
    YeastSA.sa = SA_estimated(SA_estimated~=0);
    YeastSA.grRules = model.grRules(SA_estimated~=0);
    cd SA/;
    save(['YeastSA_' expID '.mat'],'YeastSA');
    cd ../;
end

%% combine replicates
load('YeastSA_Yu_ClimR1.mat');
r1 = YeastSA;
load('YeastSA_Yu_ClimR2.mat');
r2 = YeastSA;
rxn = [r1.rxn;r2.rxn];
sa = [r1.sa;r2.sa];
grRules = [r1.grRules;r2.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_Yu_Clim.mat','YeastSA');
cd ../;

load('YeastSA_Yu_CN30R1.mat');
r1 = YeastSA;
load('YeastSA_Yu_CN30R2.mat');
r2 = YeastSA;
rxn = [r1.rxn;r2.rxn];
sa = [r1.sa;r2.sa];
grRules = [r1.grRules;r2.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_Yu_CN30.mat','YeastSA');
cd ../;

load('YeastSA_Yu_CN50R1.mat');
r1 = YeastSA;
load('YeastSA_Yu_CN50R2.mat');
r2 = YeastSA;
rxn = [r1.rxn;r2.rxn];
sa = [r1.sa;r2.sa];
grRules = [r1.grRules;r2.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_Yu_CN50.mat','YeastSA');
cd ../;

load('YeastSA_Yu_CN115R1.mat');
r1 = YeastSA;
load('YeastSA_Yu_CN115R2.mat');
r2 = YeastSA;
rxn = [r1.rxn;r2.rxn];
sa = [r1.sa;r2.sa];
grRules = [r1.grRules;r2.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_Yu_CN115.mat','YeastSA');
cd ../;

load('YeastSA_DiBartolomeo_GlucR1.mat');
r1 = YeastSA;
load('YeastSA_DiBartolomeo_GlucR2.mat');
r2 = YeastSA;
load('YeastSA_DiBartolomeo_GlucR3.mat');
r3 = YeastSA;
rxn = [r1.rxn;r2.rxn;r3.rxn];
sa = [r1.sa;r2.sa;r3.sa];
grRules = [r1.grRules;r2.grRules;r3.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_DiBartolomeo_Gluc.mat','YeastSA');
cd ../;


load('YeastSA_DiBartolomeo_EtohR1.mat');
r1 = YeastSA;
load('YeastSA_DiBartolomeo_EtohR2.mat');
r2 = YeastSA;
load('YeastSA_DiBartolomeo_EtohR3.mat');
r3 = YeastSA;
rxn = [r1.rxn;r2.rxn;r3.rxn];
sa = [r1.sa;r2.sa;r3.sa];
grRules = [r1.grRules;r2.grRules;r3.grRules];
YeastSA.rxn = unique(rxn);
YeastSA.sa = zeros(length(YeastSA.rxn),1);
YeastSA.grRules = cell(length(YeastSA.rxn),1);
for i = 1:length(YeastSA.rxn)
    YeastSA.sa(i) = mean(sa(ismember(rxn,YeastSA.rxn(i))));
    YeastSA.grRules(i) = unique(grRules(ismember(rxn,YeastSA.rxn(i))));
end
cd SA/;
save('YeastSA_DiBartolomeo_Etoh.mat','YeastSA');
cd ../;
