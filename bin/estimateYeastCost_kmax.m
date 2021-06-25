% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% cd Models/;
% save('GEM-yeast-split.mat','model_split');
% cd ../;

load('GEM-yeast-split.mat');

enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);
model = convertModel(model_split,enzymedata);

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

sol_ref = optimizeCbModel(model,'max');
gluc_ref = sol_ref.f;
model_ref = changeRxnBounds(model,substrateExID,gluc_ref,'b');
model_ref = changeObjective(model_ref, 'prot_cost_exchange');
sol_ref_min = optimizeCbModel(model_ref,'min');
flux_ref = sol_ref_min.x;

delta_gluc = zeros(1,length(AA));
delta_flux = zeros(length(model.rxns),length(AA));
delta_prot = zeros(1,length(AA));
for i = 1:length(AA)
    display(num2str(i));
    model_tmp = model;
    
    model_tmp.S(ismember(model_tmp.mets,aaList{i}),ismember(model_tmp.rxns,'r_4047')) = ...
        model_tmp.S(ismember(model_tmp.mets,aaList{i}),ismember(model_tmp.rxns,'r_4047')) - step;
    
    sol_tmp_tmp = optimizeCbModel(model_tmp,'max');
    model_tmp = changeRxnBounds(model_tmp,substrateExID,sol_tmp_tmp.f,'b');
    model_tmp = changeObjective(model_tmp, 'prot_cost_exchange');
    sol_tmp = optimizeCbModel(model_tmp,'min');
    
    idxgluc = ismember(model_tmp.rxns,substrateExID);
    
    delta_gluc(1,i) = (gluc_ref-sol_tmp_tmp.f)/mu/step;
    delta_flux(:,i) = (sol_tmp.x-flux_ref)/mu/step;
    delta_prot(1,i) = (sol_tmp.f-sol_ref_min.f)/mu/step;
end

cost_gluc = transpose(delta_gluc);
cost_prot_same = transpose(delta_prot);


cost_prot = zeros(length(AA),1);
highCSPCfraction = zeros(length(AA),1);

delta_flux(abs(delta_flux)<1e-6) = 0;
for i = 1:length(AA)
    
    flux_tmp = delta_flux(:,i);
    active_rxns = model.rxns(flux_tmp ~= 0);
    active_fluxes = flux_tmp(flux_tmp ~= 0);
    protein_cost = zeros(length(active_rxns),1);
    protein_cost_cs = zeros(length(active_rxns),1);
    
    for j = 1:length(active_rxns)
        if ismember(active_rxns(j),enzymedata.rxn)
            protein_cost(j) = enzymedata.minPC(ismember(enzymedata.rxn,active_rxns(j)));
            protein_cost_cs(j) = enzymedata.kcat_conf(ismember(enzymedata.rxn,active_rxns(j)));
        end
    end
    cost_prot(i) = sum(active_fluxes.*protein_cost)/1000;
    
    AbsPC5 = sum(abs(protein_cost(protein_cost_cs>=5).*active_fluxes(protein_cost_cs>=5)));
    totAbsPC = sum(abs(active_fluxes.*protein_cost));
    highCSPCfraction(i) = AbsPC5 / totAbsPC;
end
cost_kmax.AA = AA;
cost_kmax.cost_gluc = cost_gluc;
cost_kmax.cost_prot = cost_prot;
cost_kmax.highCSPCfraction = highCSPCfraction;
cd YeastCost/;
save('cost_kmax.mat','cost_kmax');
cd ../;