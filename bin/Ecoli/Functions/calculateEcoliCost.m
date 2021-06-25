function cost_ecoli = calculateEcoliCost(model,substrateExID,mu,step,tol)

aaList = {'ala__L[c]';'arg__L[c]';'asn__L[c]';'asp__L[c]';'cys__L[c]';
'gln__L[c]';'glu__L[c]';'gly[c]';'his__L[c]';'ile__L[c]';
'leu__L[c]';'lys__L[c]';'met__L[c]';'phe__L[c]';'pro__L[c]';
'ser__L[c]';'thr__L[c]';'trp__L[c]';'tyr__L[c]';'val__L[c]'};
AA = {'A';'R';'N';'D';'C';
    'Q';'E';'G';'H';'I';
    'L';'K';'M';'F';'P';
    'S';'T';'W';'Y';'V'};

model = changeObjective(model, substrateExID);
model = changeRxnBounds(model, substrateExID, -1000, 'l'); %glucose exchange
model = changeRxnBounds(model, substrateExID, 0, 'u'); %glucose exchange
model = changeRxnBounds(model, 'EX_co2_e', 0, 'l');
model = changeRxnBounds(model, 'prot_cost_exchange', 0, 'l');
model = changeRxnBounds(model, 'prot_cost_exchange', 100000, 'u');
model = changeRxnBounds(model, 'ATPM', 6.86, 'b');
model = changeRxnBounds(model, 'BIOMASS_Ec_iML1515_core_75p37M', mu, 'b'); %growth rate

% change biomass to protein
atpBiomass = {'h[c]';'adp[c]';'pi[c]';'atp[c]';'h2o[c]'};
model.S(ismember(model.mets,atpBiomass),ismember(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M')) = ...
    [75.37;75.37;75.37;-75.37;-75.37];
keyBiomass = [aaList;atpBiomass];
allBiomass = model.mets(find(model.S(:,ismember(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M'))));
otherBiomass = allBiomass(~ismember(allBiomass,keyBiomass));
model.S(ismember(model.mets,otherBiomass),ismember(model.rxns,'BIOMASS_Ec_iML1515_core_75p37M')) = 0;

sol_org_1 = optimizeCbModel(model,'max','one');
model_org_1 = changeRxnBounds(model,substrateExID,sol_org_1.f-tol,'b');
model_org_1 = changeObjective(model_org_1, 'prot_cost_exchange');
sol_org_new_1 = optimizeCbModel(model_org_1,'min');
tot_prot_cost_org_1 = sol_org_new_1.f;

flux_ref_1 = sol_org_new_1.x;

cost_gluc_1 = zeros(length(AA),1);
cost_prot_1 = zeros(length(AA),1);
flux_AAs_1 = zeros(length(model.rxns),length(AA));

for i = 1:length(aaList)
    display(num2str(i));
    model_tmp_1 = model;
    model_tmp_1.S(ismember(model_tmp_1.mets,aaList{i}),ismember(model_tmp_1.rxns,'BIOMASS_Ec_iML1515_core_75p37M')) = ...
        model_tmp_1.S(ismember(model_tmp_1.mets,aaList{i}),ismember(model_tmp_1.rxns,'BIOMASS_Ec_iML1515_core_75p37M')) - step;
    sol_tmp_1 = optimizeCbModel(model_tmp_1,'max','one');
    
    model_tmp_tmp_1 = changeRxnBounds(model_tmp_1,substrateExID,sol_tmp_1.f-tol,'b');
    model_tmp_tmp_1 = changeObjective(model_tmp_tmp_1, 'prot_cost_exchange');
    sol_tmp_tmp_1 = optimizeCbModel(model_tmp_tmp_1,'min');
    
    flux_AAs_1(:,i) = sol_tmp_tmp_1.x;
    
    idxgluc = ismember(model.rxns,substrateExID);
    cost_gluc_1(i) = (sol_org_new_1.x(idxgluc)-sol_tmp_tmp_1.x(idxgluc))/step/mu;
    cost_prot_1(i) = (sol_tmp_tmp_1.f-tot_prot_cost_org_1)/step/mu;
    
end

cost_ecoli.AA = AA;
cost_ecoli.cost_gluc = cost_gluc_1;
cost_ecoli.cost_prot = cost_prot_1;
cost_ecoli.flux_ref = flux_ref_1;
cost_ecoli.flux_AAs = flux_AAs_1;






