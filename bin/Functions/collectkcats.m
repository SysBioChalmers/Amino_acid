%% collectkcats
function enzymedata = collectkcats(rxnlist,model,org_name)
%% Load kcats data from BRENDA database
allkcats = loadkcats;

%% Collect EC number
ecexcel = strcat('ec_',strrep(org_name,' ','_'),'.xlsx');
[~, rawUniprot, ~] = xlsread(ecexcel,'Uniprot');
[~, rawKEGG, ~] = xlsread(ecexcel,'KEGG');

id_Uniprot = rawUniprot(:,1);
ec_Uniprot = rawUniprot(:,2);
id_KEGG = rawKEGG(:,1);
ec_KEGG = rawKEGG(:,2);

id_list = unique([id_Uniprot;id_KEGG]);

ecdata = struct();
ecdata.id = cell(0,1);
ecdata.ec = cell(0,1);

for i = 1:length(id_list)
    id = id_list(i);
    if ismember(id,id_Uniprot)
        ec_U_tmp = ec_Uniprot(ismember(id_Uniprot,id));
        ec_U_tmp = split(ec_U_tmp);
        id_U_tmp = repelem(id,length(ec_U_tmp))';
        ecdata.id = [ecdata.id;id_U_tmp];
        ecdata.ec = [ecdata.ec;ec_U_tmp];
    end
    if ismember(id,id_KEGG)
        ec_K_tmp = ec_KEGG(ismember(id_KEGG,id));
        id_K_tmp = repelem(id,length(ec_K_tmp))';
        ecdata.id = [ecdata.id;id_K_tmp];
        ecdata.ec = [ecdata.ec;ec_K_tmp];
    end
end

%% Assign kcat for each reaction with GPR

enzymedata = struct();
enzymedata.rxn = rxnlist;
enzymedata.substrate = cell(length(rxnlist),1);
enzymedata.subunit = cell(length(rxnlist),25);
enzymedata.subunit_ec = cell(length(rxnlist),25);
enzymedata.subunit_kcat = zeros(length(rxnlist),25);
enzymedata.subunit_kcat_conf = nan(length(rxnlist),25);
enzymedata.kcat = zeros(length(rxnlist),1);
enzymedata.kcat_conf = zeros(length(rxnlist),1);

for i = 1:length(enzymedata.rxn)
    
    disp(['Assigning kcats:' num2str(i) '/' num2str(length(enzymedata.rxn))]);
    
    rxn_id = enzymedata.rxn{i};
    
    % add substrates
    idx_tmp = ismember(model.rxns,rxn_id);
    s_tmp = model.S(:,idx_tmp);
    mets_tmp = model.metNames(s_tmp < 0);
    if strcmp(org_name,'saccharomyces cerevisiae')
        mets_tmp = cellfun(@(x) x(1:strfind(x,' [')-1),mets_tmp,'UniformOutput',false);
    end
    mets_tmp = strjoin(mets_tmp,'; ');
    enzymedata.substrate(i,1) = {mets_tmp};
    
    % add subunits
    gr_tmp = model.grRules{ismember(model.rxns,rxn_id)};
    gr_tmp = strrep(gr_tmp,'(','');
    gr_tmp = strrep(gr_tmp,')','');
    gr_tmp = strrep(gr_tmp,' and ',' ');
    gr_tmp = strrep(gr_tmp,' or ',' ');
    gr_tmp = strtrim(gr_tmp);
    gr_tmp = strsplit(gr_tmp);
    gr_tmp = unique(gr_tmp);
    na_tmp = repelem({''},25-length(gr_tmp));
    enzymedata.subunit(i,:) = [gr_tmp na_tmp];
    
    % add ec numbers of subunits
    for j = 1:length(gr_tmp)
        subunit_tmp = gr_tmp(j);
        if ismember(subunit_tmp,ecdata.id)
            ec_tmp = unique(ecdata.ec(ismember(ecdata.id,subunit_tmp)));
            ec_tmp = strjoin(ec_tmp,'; ');
        else
            ec_tmp = 'NA';
        end
        enzymedata.subunit_ec(i,j) = {ec_tmp};
    end
    
    % assign kcats and confidence scores for subunits
    for j = 1:length(gr_tmp)
        ec_tmp = enzymedata.subunit_ec(i,j);
        substrate_tmp = enzymedata.substrate(i,1);
        [finalkcat_tmp, conf_score_tmp] = searchkcat(ec_tmp,substrate_tmp,org_name,allkcats);
        enzymedata.subunit_kcat(i,j) = finalkcat_tmp;
        enzymedata.subunit_kcat_conf(i,j) = conf_score_tmp;
    end
    
    % assign kcats and confidence scores for enzymes
    conf_tmp = max(enzymedata.subunit_kcat_conf(i,:));
    enzymedata.kcat_conf(i,1) = conf_tmp;
    finalkcat_list = enzymedata.subunit_kcat(i,:);
    kcat_tmp = finalkcat_list(enzymedata.subunit_kcat_conf(i,:) == conf_tmp);
    enzymedata.kcat(i,1) = median(kcat_tmp); %choose median among subunits/isozymes
end

% no kcat assigned enzyme assumed to be the median of the collected dataset
medianvalue = median(enzymedata.kcat(~isnan(enzymedata.kcat)));
enzymedata.kcat(isnan(enzymedata.kcat)) = medianvalue;



    