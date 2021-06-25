%% splitModel 
function new_model = splitRevRxns(model)

disp('Splitting reactions...');
% Change gene ID in the model.rules into that in the model.genes.
for i = 1:length(model.genes)
    old = strcat('x(',num2str(i),')');
    new = model.genes{i};
    model.rules = cellfun(@(x) strrep(x,old,new),...
                            model.rules,'UniformOutput',false);
end

% Generate a matrix for splitting reversible reactions.
matrix = struct();
matrix.RxnList = cell(0,1);
matrix.CompoList = cell(0,1);
matrix.CoeffList = zeros(0,1);
matrix.CatalystList = cell(0,1);
matrix.LBList = zeros(0,1);
matrix.UBList = zeros(0,1);

UnqRxnList = model.rxns;
for i = 1:length(UnqRxnList)
%     id_tmp = UnqRxnList(i);
%     idx = ismember(matrix1.RxnList,id_tmp);
    
    idx_substrate = model.S(:,i) < 0;
    idx_product = model.S(:,i) > 0;
    CompoS = model.mets(idx_substrate);
    CompoP = model.mets(idx_product);
    CoeffS = model.S(idx_substrate,i);
    CoeffP = model.S(idx_product,i);
    CompoList = [CompoS;CompoP];
    CoeffList = [CoeffS;CoeffP];
    n = length(CompoList);
    Rxn = model.rxns(i);
    Catalyst = model.rules(i);
    lb_tmp = model.lb(i);
    ub_tmp = model.ub(i);
    RxnList = repmat(Rxn,n,1);
    CatalystList = repmat(Catalyst,n,1);
    LBList = repmat(lb_tmp,n,1);
    UBList = repmat(ub_tmp,n,1);
    
%     RxnList = matrix1.RxnList(idx);
%     CompoList = matrix1.CompoList(idx);
%     CoeffList = matrix1.CoeffList(idx);
%     CatalystList = matrix1.CatalystList(idx);
%     LBList = matrix1.LBList(idx);
%     UBList = matrix1.UBList(idx);
%     n = length(RxnList);
    
    x = length(matrix.RxnList);%count rows in matrix2
    if isempty(CatalystList{1})%if is spontaneous reaction
        matrix.RxnList(x+1:x+n,1) = RxnList;
        matrix.CompoList(x+1:x+n,1) = CompoList;
        matrix.CoeffList(x+1:x+n,1) = CoeffList;
        matrix.CatalystList(x+1:x+n,1) = CatalystList;
        matrix.LBList(x+1:x+n,1) = LBList;
        matrix.UBList(x+1:x+n,1) = UBList;
    else%if is enzymatic reaction
        if LBList(1) >= 0%if is not reversible
            matrix.RxnList(x+1:x+n,1) = RxnList;
            matrix.CompoList(x+1:x+n,1) = CompoList;
            matrix.CoeffList(x+1:x+n,1) = CoeffList;
            matrix.CatalystList(x+1:x+n,1) = CatalystList;
            matrix.LBList(x+1:x+n,1) = LBList;
            matrix.UBList(x+1:x+n,1) = UBList;
        else%if is reversible
            %add forward rows
            if UBList(1) > 0
                rxnname_tmp = strcat(RxnList{1},'_fwd');
                matrix.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
                matrix.CompoList(x+1:x+n,1) = CompoList;
                matrix.CoeffList(x+1:x+n,1) = CoeffList;
                matrix.CatalystList(x+1:x+n,1) = CatalystList;
                matrix.LBList(x+1:x+n,1) = zeros(n,1);
                matrix.UBList(x+1:x+n,1) = UBList;
            end
            x = length(matrix.RxnList);
            %add reverse rows
            rxnname_tmp = strcat(RxnList{1},'_rvs');
            matrix.RxnList(x+1:x+n,1) = repmat({rxnname_tmp},n,1);
            matrix.CompoList(x+1:x+n,1) = CompoList;
            matrix.CoeffList(x+1:x+n,1) = -1*CoeffList;
            matrix.CatalystList(x+1:x+n,1) = CatalystList;
            matrix.LBList(x+1:x+n,1) = zeros(n,1);
            matrix.UBList(x+1:x+n,1) = -1*LBList;
        end
    end
end

matrixSplit = matrix;

% Converted to COBRA model
new_model = struct();
new_model.rxns = cell(0,1);
new_model.lb = zeros(0,1);
new_model.ub = zeros(0,1);
new_model.mets = cell(0,1);
new_model.S = sparse(0,0);
new_model.b = zeros(0,1);
new_model.rxnGeneMat = sparse(0,0);
new_model.c = zeros(0,1);
new_model.rules = cell(0,1);
new_model.genes = cell(0,1);
new_model.csense = char();

new_model.metFormulas = cell(0,1);
new_model.metNames = cell(0,1);
new_model.rxnNames = cell(0,1);

new_model.osenseStr = model.osenseStr;

if isfield(model,'compNames')
    new_model.compNames = model.compNames;
end
if isfield(model,'comps')
    new_model.comps = model.comps;
end

UnqRxnList = unique(matrixSplit.RxnList);
for i = 1:length(UnqRxnList)
    
    id_tmp = UnqRxnList(i);
    idx = ismember(matrixSplit.RxnList,id_tmp);
    
    RxnList = matrixSplit.RxnList(idx);
    CompoList = matrixSplit.CompoList(idx);
    CoeffList = matrixSplit.CoeffList(idx);
    CatalystList = matrixSplit.CatalystList(idx);
    LBList = matrixSplit.LBList(idx);
    UBList = matrixSplit.UBList(idx);
    
    rxn_tmp = RxnList{1};
    catalyst_tmp = CatalystList{1};
    lb_tmp = unique(LBList);
    ub_tmp = unique(UBList);
    
    % add metabolites
    [isInModelidx,~] = ismember(CompoList,new_model.mets);
    NotInModel = CompoList(~isInModelidx);
    idx_met_tmp = ismember(model.mets,NotInModel);
	metid_tmp = model.mets(idx_met_tmp);
	metname_tmp = model.metNames(idx_met_tmp);
	metformula_tmp = model.metFormulas(idx_met_tmp);
    new_model = addMetabolite(new_model,metid_tmp,...
                              'metName',metname_tmp,...
                              'metFormula',metformula_tmp);
    
    % add reactions
    if length(rxn_tmp) > 5
        if strcmp('_fwd',rxn_tmp(end-3:end)) || strcmp('_rvs',rxn_tmp(end-3:end))
            newrxntmp = rxn_tmp(1:end-4);
        else
            newrxntmp = rxn_tmp;
        end
    else
        newrxntmp = rxn_tmp;
    end
    
    reactionname_tmp = model.rxnNames{ismember(model.rxns,newrxntmp)};
    new_model = addReaction(new_model,rxn_tmp,...
                            'reactionName',reactionname_tmp,...
                            'metaboliteList',CompoList,...
                            'stoichCoeffList',CoeffList,...
                            'lowerBound',lb_tmp,...
                            'upperBound',ub_tmp,...
                            'geneRule',catalyst_tmp);
end
end

