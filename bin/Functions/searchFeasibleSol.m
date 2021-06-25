function [fluxes, deviation] = searchFeasibleSol(model,exRxnList,exFluxes,osenseStr,step)
% This is for the case that no solution can be obtained for a given model
% with experimentally measured data as constraint. In order to get a
% feasible solution, we adjust exchange reaction rates gradually.

model_org = changeRxnBounds(model,exRxnList,exFluxes,'b');
sol = optimizeCbModel(model_org,osenseStr,'one');

deviation = step;
while ~strcmp(sol.origStat,'OPTIMAL')
    exFluxes_lb = exFluxes;
    exFluxes_ub = exFluxes;
    
    exFluxes_lb(exFluxes_lb<0) = exFluxes_lb(exFluxes_lb<0)*(1+deviation);
    exFluxes_lb(exFluxes_lb>0) = exFluxes_lb(exFluxes_lb>0)*(1-deviation);
    
    exFluxes_ub(exFluxes_ub<0) = exFluxes_ub(exFluxes_ub<0)*(1-deviation);
    exFluxes_ub(exFluxes_ub>0) = exFluxes_ub(exFluxes_ub>0)*(1+deviation);
    
    model_tmp = model;
    model_tmp = changeRxnBounds(model_tmp,exRxnList,exFluxes_lb,'l');
    model_tmp = changeRxnBounds(model_tmp,exRxnList,exFluxes_ub,'u');
    sol = optimizeCbModel(model_tmp,osenseStr,'one');
    
    deviation = deviation + step;
end

fluxes = sol.x;