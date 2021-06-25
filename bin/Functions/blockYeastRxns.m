%% blockYeastRxns 
function model = blockYeastRxns(model)
% Some reactions should be block to avoid weird flux distributions.

% block Gcy1, an alternative glycerol dissimilation pathway that is active 
% under microaerobic conditions (PMID: 22979944)
model = changeRxnBounds(model,'r_0487',0,'b');

model = changeRxnBounds(model,'r_4262_fwd',0,'b'); % citrate hydroxymutase
model = changeRxnBounds(model,'r_4262_rvs',0,'b'); % citrate hydroxymutase

% block some reactions that done in PMID: 28779005.
model = changeRxnBounds(model,'r_2045_rvs',0,'b'); % serine transport from [m] to [c]
model = changeRxnBounds(model,'r_0659_fwd',0,'b'); % isocitrate dehydrogenase (NADP)
model = changeRxnBounds(model,'r_0659_rvs',0,'b'); % isocitrate dehydrogenase (NADP)

% model = changeRxnBounds(model,'r_0725_fwd',0,'b'); % methenyltetrahydrofolate cyclohydrolase
% model = changeRxnBounds(model,'r_0918',0,'b'); % phosphoserine transaminase

model = changeRxnBounds(model,'r_4216_rvs',0,'b'); % block the reaction to produce FMN without ATP

% The following two reactions account for the transport of pyruvate from
% [c] to [m] without enzyme cost, should be blocked.
model = changeRxnBounds(model,'r_1137',0,'b');
model = changeRxnBounds(model,'r_1138',0,'b');

% Block backward reactions for PLP synthesis.
plpbwrxns = {'r_4211_rvs' 'r_4212_rvs'};
model = changeRxnBounds(model,plpbwrxns,zeros(1,length(plpbwrxns)),'b');


