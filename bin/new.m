% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% cd Models/;
% save('GEM-yeast-split.mat','model_split');
% cd ../;

load('GEM-yeast-split.mat');

% enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
% enzymedata = updateYeastEnzyme(enzymedata,model_split);
% model = convertModel(model_split,enzymedata);

model = model_split;
model = setYeastMedia(model);
% model = blockYeastRxns(model);%%%%%%%%%%%%%%


step = 0.001;
mu = 1;
substrateExID = 'r_1714'; % glucose

aasubsList = {'s_0404[c]';'s_0428[c]';'s_0430[c]';'s_0432[c]';'s_0542[c]';
's_0747[c]';'s_0748[c]';'s_0757[c]';'s_0832[c]';'s_0847[c]';
's_1077[c]';'s_1099[c]';'s_1148[c]';'s_1314[c]';'s_1379[c]';
's_1428[c]';'s_1491[c]';'s_1527[c]';'s_1533[c]';'s_1561[c]'};
aaprodList = {'s_1582[c]';'s_1583[c]';'s_1585[c]';'s_1587[c]';'s_1589[c]';
's_1590[c]';'s_1591[c]';'s_1593[c]';'s_1594[c]';'s_1596[c]';
's_1598[c]';'s_1600[c]';'s_1602[c]';'s_1604[c]';'s_1606[c]';
's_1607[c]';'s_1608[c]';'s_1610[c]';'s_1612[c]';'s_1614[c]'};
aaList = {'s_0955[c]';'s_0965[c]';'s_0969[c]';'s_0973[c]';'s_0981[c]';
's_0999[c]';'s_0991[c]';'s_1003[c]';'s_1006[c]';'s_1016[c]';
's_1021[c]';'s_1025[c]';'s_1029[c]';'s_1032[c]';'s_1035[c]';
's_1039[c]';'s_1045[c]';'s_1048[c]';'s_1051[c]';'s_1056[c]'};

AA = {'A';'R';'N';'D';'C';
    'Q';'E';'G';'H';'I';
    'L';'K';'M';'F';'P';
    'S';'T';'W';'Y';'V'};

% change biomass to protein
otherBiomass = {'s_1096[c]';'s_3718[c]';'s_3719[c]';'s_3720[c]';'s_4205[c]';'s_4206[c]'};
model.S(ismember(model.mets,otherBiomass),ismember(model.rxns,'r_4041')) = 0;

% remove tRNA
[~,p_aasubs] = ismember(aasubsList,model.mets);
aacoef = model.S(p_aasubs,ismember(model.rxns,'r_4047'));
[~,p_aa] = ismember(aaList,model.mets);
[~,p_aaprod] = ismember(aaprodList,model.mets);
model.S(p_aasubs,ismember(model.rxns,'r_4047')) = 0;
model.S(p_aa,ismember(model.rxns,'r_4047')) = aacoef;
model.S(p_aaprod,ismember(model.rxns,'r_4047')) = 0;

model = changeObjective(model, substrateExID);
model = changeRxnBounds(model, substrateExID, -1000, 'l');
model = changeRxnBounds(model, substrateExID, 0, 'u');

model = changeRxnBounds(model, 'r_2111', mu, 'b'); %growth rate

sol_ref = optimizeCbModel(model,'max','one');
gluc_ref = sol_ref.f;
flux_ref = sol_ref.x;

delta_gluc = zeros(1,length(AA));
delta_flux = zeros(length(model.rxns),length(AA));
for i = 1:length(AA)
    display(num2str(i));
    model_tmp = model;
    
    model_tmp.S(ismember(model_tmp.mets,aaList{i}),ismember(model_tmp.rxns,'r_4047')) = ...
        model_tmp.S(ismember(model_tmp.mets,aaList{i}),ismember(model_tmp.rxns,'r_4047')) - step;
    
    sol_tmp = optimizeCbModel(model_tmp,'max','one');
    
%     sol_tmp_tmp = optimizeCbModel(model_tmp,'max','one');
%     model_tmp = changeRxnBounds(model_tmp,substrateExID,sol_tmp_tmp.f,'b');
%     [sol_tmp, ~, ~, ~] = MOMA(model, model_tmp, 'max', false, 'one');
    
    gluc_tmp = sol_tmp.f;
    delta_gluc(1,i) = abs((sol_tmp.f-gluc_ref))/mu/step;
    delta_flux(:,i) = (sol_tmp.x-flux_ref)/mu/step;
end

cost_gluc = transpose(delta_gluc);
cost_prot = zeros(length(AA),1);

enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);
delta_flux(abs(delta_flux)<1e-6) = 0;
for i = 1:length(AA)
    
    flux_tmp = delta_flux(:,i);
    active_rxns = model.rxns(flux_tmp ~= 0);
    active_fluxes = flux_tmp(flux_tmp ~= 0);
    
    tot_prot_cost = 0;
    for j = 1:length(active_rxns)
        if ismember(active_rxns(j),enzymedata.rxn)
            prot_cost_rxn = enzymedata.minPC(ismember(enzymedata.rxn,active_rxns(j)))*active_fluxes(j);
            tot_prot_cost = tot_prot_cost + prot_cost_rxn;
        end
    end
    cost_prot(i) = tot_prot_cost/1000;
end
cost_maxSA.AA = AA;
cost_maxSA.cost_gluc = cost_gluc;
cost_maxSA.cost_prot = cost_prot;
cd YeastCost/;
save('cost_maxSA.mat','cost_maxSA');
cd ../;

% condition-specific sa for estimating cost
load('YeastSA.mat');
cond = YeastSA.condition(~contains(YeastSA.condition,'etoh','IgnoreCase',true));

for k = 1:length(cond)
    display([num2str(k),'/',num2str(length(cond))]);
    enzymedata_cond = enzymedata;
    enzymedata_cond = updateCondSpecYeastEnzyme(enzymedata_cond,cond(k));
    
    for i = 1:length(AA)
        
        flux_tmp = delta_flux(:,i);
        active_rxns = model.rxns(flux_tmp ~= 0);
        active_fluxes = flux_tmp(flux_tmp ~= 0);
        
        tot_prot_cost = 0;
        for j = 1:length(active_rxns)
            if ismember(active_rxns(j),enzymedata_cond.rxn)
                prot_cost_rxn = enzymedata_cond.minPC(ismember(enzymedata_cond.rxn,active_rxns(j)))*active_fluxes(j);
                tot_prot_cost = tot_prot_cost + prot_cost_rxn;
            end
        end
        cost_prot(i) = tot_prot_cost/1000;
    end
    cost_CondSA.AA = AA;
    cost_CondSA.cost_gluc = cost_gluc;
    cost_CondSA.cost_prot = cost_prot;
    cd YeastCost/;
    save(['cost_CondSA_' cond{k} '.mat'],'cost_CondSA');
    cd ../;
end



