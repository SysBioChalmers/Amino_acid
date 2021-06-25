%% costPerProtein 
function [cost_total, cost_average] = costPerProtein(gene_list,AA,AA_cost,org_name)
% calculate cost for each protein in a given list

[~,txt,~] = xlsread(strcat('UniProt_',org_name,'.xlsx'));
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);

cost_total = zeros(length(gene_list),1);
seq_length = zeros(length(gene_list),1);

for i = 1:length(gene_list)
    
    if ismember(gene_list(i),seq_data.gene)
        seq = seq_data.seq{ismember(seq_data.gene,gene_list(i))};
        seq_length(i,1) = length(seq);
        AAcount = zeros(length(AA),1);
        for j = 1:length(AA)
            AAcount(j,1) = length(strfind(seq,AA{j}));
        end
        cost_total(i,1) = sum(AAcount.*AA_cost);
    else
        seq_length(i,1) = 0;
        cost_total(i,1) = 0;
    end
end

cost_average = cost_total./seq_length;

cost_average(isnan(cost_average)) = 0;