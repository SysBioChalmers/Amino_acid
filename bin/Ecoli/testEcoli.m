clear;
load('GEM-ecoli-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'escherichia coli');
enzymedata = updateEcoliEnzyme(enzymedata,model_split);

load('EcoliSA.mat');

% k = 28;
% enzymedata = updateCondSpecEcoliEnzyme(enzymedata,EcoliSA.condition(k));

model = convertModel(model_split,enzymedata);
model = blockEcoliRxns(model);



model = changeRxnBounds(model, 'EX_glc__D_e', -1, 'b'); % gluc
model = changeRxnBounds(model, 'ATPM', 0, 'l'); % ATPM
model = changeRxnBounds(model, 'ATPM', 1000, 'u'); % ATPM

model = changeRxnBounds(model, 'prot_cost_exchange', 0, 'l');
model = changeRxnBounds(model, 'prot_cost_exchange', 1000, 'u');

modelly = changeObjective(model, 'ATPM'); % ac
modelly = changeRxnBounds(modelly, 'EX_ac_e', 2, 'b'); % ac
solly = optimizeCbModel(modelly,'max','one');
pcly = solly.x(ismember(modelly.rxns,'prot_cost_exchange'));
pely = solly.x(ismember(modelly.rxns,'ATPM'))/pcly;
rxnly = enzymedata.rxn(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));
rxnpcly = enzymedata.minPC(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));
rxnconfly = enzymedata.kcat_conf(ismember(enzymedata.rxn,modelly.rxns(solly.x ~= 0)));

modelhy = changeObjective(model, 'ATPM'); % ATP
solhy = optimizeCbModel(modelhy,'max','one');
pchy = solhy.x(ismember(modelhy.rxns,'prot_cost_exchange'));
pehy = solhy.x(ismember(modelhy.rxns,'ATPM'))/pchy;
rxnhy = enzymedata.rxn(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));
rxnpchy = enzymedata.minPC(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));
rxnconfhy = enzymedata.kcat_conf(ismember(enzymedata.rxn,modelhy.rxns(solhy.x ~= 0)));

ratio = pely/pehy
