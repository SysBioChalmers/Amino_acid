load('GEM-yeast.mat');
for i = 1:length(model.genes)
    old = strcat('x(',num2str(i),')');
    new = model.genes{i};
    model.rules = cellfun(@(x) strrep(x,old,new),...
                            model.rules,'UniformOutput',false);
end
clear i new old;
idx = contains(model.rules,' | ') & ~contains(model.rules,' & ');
iso_enzymes = model.rules(idx);
iso_rxns = model.rxns(idx);
unq_isozyme_groups = unique(iso_enzymes);

load('cost_kmax.mat')
[num,txt,~] = xlsread('paxdb_yeast.xlsx'); 
abund = num;
protID = txt(2:end,1);
clear num txt;

res_tf = false(0,2);
res_rxnlist = cell(0,1);
for i = 1:length(unq_isozyme_groups)
    
    disp([num2str(i),'/',num2str(length(unq_isozyme_groups))]);
    
    tmp = unq_isozyme_groups(i);
    tmp = strrep(tmp,'( ','');
    tmp = strrep(tmp,' )','');
    isozymes = split(tmp,' | ');
    
    ovlp_isozymes = intersect(isozymes,protID);
    if length(ovlp_isozymes) > 1
        [~,b] = ismember(ovlp_isozymes,protID);
        abund_isozymes = abund(b);

        gluc_cost_isozymes = calculateCost4Protein(ovlp_isozymes,cost_kmax.AA,cost_kmax.cost_gluc,'Yeast',1);
        prot_cost_isozymes = calculateCost4Protein(ovlp_isozymes,cost_kmax.AA,cost_kmax.cost_prot,'Yeast',1);

        max_abd_idx = abund_isozymes == max(abund_isozymes);
        g_cost_max = mean(gluc_cost_isozymes(max_abd_idx));
        p_cost_max = mean(prot_cost_isozymes(max_abd_idx));

        min_abd_idx = abund_isozymes == min(abund_isozymes);
        g_cost_min = mean(gluc_cost_isozymes(min_abd_idx));
        p_cost_min = mean(prot_cost_isozymes(min_abd_idx));

        res_tf = [res_tf;[g_cost_max < g_cost_min,p_cost_max < p_cost_min]];
        
        res_rxnlist = [res_rxnlist;join(iso_rxns(ismember(iso_enzymes,unq_isozyme_groups(i))),';')];
    end
end


% 
% 
% 
% ranksump = ranksum([1.005 1.01],[1.2 1.16])
% 
% [~,p1] = ttest2([3.2 2.6 3.5],[0.2 0.6 0.5],'Vartype','unequal')

