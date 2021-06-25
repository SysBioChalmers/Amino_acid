clear;
load('GEM-yeast-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);

load('YeastSA.mat');

% k = 4;
% % load(['../Yeast_kapp/kappEstimation/Fluxes/Fluxes_' YeastSA.condition{k} '.mat']);
% % zerorxns = Fluxes.model.rxns(Fluxes.pFBA == 0);
% enzymedata = updateCondSpecYeastEnzyme(enzymedata,YeastSA.condition(k));

model = convertModel(model_split,enzymedata);
model = setYeastMedia(model);
model = blockYeastRxns(model);

% model.lb(ismember(model.rxns,zerorxns)) = 0;
% model.ub(ismember(model.rxns,zerorxns)) = 0;



model = changeRxnBounds(model, 'r_1714', -1, 'b'); % gluc
model = changeRxnBounds(model, 'r_4046', 0, 'l'); % ATPM
model = changeRxnBounds(model, 'r_4046', 1000, 'u'); % ATPM

model = changeRxnBounds(model, 'prot_cost_exchange', 0, 'l');
model = changeRxnBounds(model, 'prot_cost_exchange', 1000, 'u');

modelly = changeObjective(model, 'r_1761'); % etoh
solly = optimizeCbModel(modelly,'max','one');
pcly = solly.x(ismember(modelly.rxns,'prot_cost_exchange'));
pely = solly.x(ismember(modelly.rxns,'r_4046'))/pcly;
rxnly = enzymedata.rxn(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));
rxnpcly = enzymedata.minPC(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));
rxnconfly = enzymedata.kcat_conf(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));

% modelhy = changeObjective(model, 'r_1672'); % CO2
modelhy = changeObjective(model, 'r_4046'); % ATP
solhy = optimizeCbModel(modelhy,'max','one');
pchy = solhy.x(ismember(modelhy.rxns,'prot_cost_exchange'));
pehy = solhy.x(ismember(modelhy.rxns,'r_4046'))/pchy;
rxnhy = enzymedata.rxn(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));
rxnpchy = enzymedata.minPC(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));
rxnconfhy = enzymedata.kcat_conf(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));

ratio = pely/pehy
