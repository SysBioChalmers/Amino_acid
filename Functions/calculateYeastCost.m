function cost_yeast = calculateYeastCost(model,substrateExID,mu,step)

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

model = changeRxnBounds(model, 'prot_cost_exchange', 0, 'l');
model = changeRxnBounds(model, 'prot_cost_exchange', 100000, 'u');

model = changeRxnBounds(model, 'r_2111', mu, 'b'); %growth rate

sol_org_1 = optimizeCbModel(model,'max');
model_org_1 = changeRxnBounds(model,substrateExID,sol_org_1.f,'b');
model_org_1 = changeObjective(model_org_1, 'prot_cost_exchange');
sol_org_new_1 = optimizeCbModel(model_org_1,'min');
tot_prot_cost_org_1 = sol_org_new_1.f;

flux_ref_1 = sol_org_new_1.x;

cost_gluc_1 = zeros(length(AA),1);
cost_prot_1 = zeros(length(AA),1);
flux_AAs_1 = zeros(length(model.rxns),length(AA));
for i = 1:length(AA)
    display(num2str(i));
    model_tmp_1 = model;
    
    model_tmp_1.S(ismember(model_tmp_1.mets,aaList{i}),ismember(model_tmp_1.rxns,'r_4047')) = ...
        model_tmp_1.S(ismember(model_tmp_1.mets,aaList{i}),ismember(model_tmp_1.rxns,'r_4047')) - step;
    
    sol_tmp_1 = optimizeCbModel(model_tmp_1,'max');
    
    model_tmp_tmp_1 = changeRxnBounds(model_tmp_1,substrateExID,sol_tmp_1.f,'b');
    model_tmp_tmp_1 = changeObjective(model_tmp_tmp_1, 'prot_cost_exchange');
    sol_tmp_tmp_1 = optimizeCbModel(model_tmp_tmp_1,'min');
    
    if strcmp(sol_tmp_tmp_1.origStat,'OPTIMAL')
        flux_AAs_1(:,i) = sol_tmp_tmp_1.x;
        idxgluc = ismember(model.rxns,substrateExID);
        cost_gluc_1(i) = (sol_org_new_1.x(idxgluc)-sol_tmp_tmp_1.x(idxgluc))/step/mu;
        cost_prot_1(i) = (sol_tmp_tmp_1.f-tot_prot_cost_org_1)/step/mu;
    end
end

cost_yeast.AA = AA;
cost_yeast.cost_gluc = cost_gluc_1;
cost_yeast.cost_prot = cost_prot_1;
cost_yeast.flux_ref = flux_ref_1;
cost_yeast.flux_AAs = flux_AAs_1;
