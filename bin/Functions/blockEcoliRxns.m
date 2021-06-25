%% blockEcoliRxns 
function model = blockEcoliRxns(model)
% Some reactions should be block to avoid weird flux distributions.

% block reactions as in iJO1366
model = changeRxnBounds(model,'DHPTDNR',0,'b');
model = changeRxnBounds(model,'DHPTDNRN',0,'b');
model = changeRxnBounds(model,'SUCASPtpp_fwd',0,'b');
model = changeRxnBounds(model,'SUCASPtpp_rvs',0,'b');
model = changeRxnBounds(model,'SUCFUMtpp_fwd',0,'b');
model = changeRxnBounds(model,'SUCFUMtpp_rvs',0,'b');
model = changeRxnBounds(model,'SUCMALtpp_fwd',0,'b');
model = changeRxnBounds(model,'SUCMALtpp_rvs',0,'b');
model = changeRxnBounds(model,'SUCTARTtpp_fwd',0,'b');
model = changeRxnBounds(model,'SUCTARTtpp_rvs',0,'b');
model = changeRxnBounds(model,'CAT',0,'b');
model = changeRxnBounds(model,'FHL',0,'b');
model = changeRxnBounds(model,'SPODM',0,'b');
model = changeRxnBounds(model,'SPODMpp',0,'b');
% block other glc transporters
model = changeRxnBounds(model,'GLCt2pp',0,'b');
model = changeRxnBounds(model,'GLCabcpp',0,'b');
% block others
model = changeRxnBounds(model,'PFL',0,'b'); % only active under anaerobic condition
model = changeRxnBounds(model,'DHAPT',0,'b'); % growth on glycerol