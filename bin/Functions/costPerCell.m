%% costPerCell
function cost = costPerCell(gene_list_measured,abund_list_measured,AA,AA_cost,org_name)
% calculate cost for a cell

[~, cost_gene] = costPerProtein(gene_list_measured,AA,AA_cost,org_name);


[~,txt,~] = xlsread(strcat('UniProt_',org_name,'.xlsx'));
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);

length_gene = zeros(length(gene_list_measured),1);
for i = 1:length(gene_list_measured)
    if ismember(gene_list_measured(i),seq_data.gene)
        seq = seq_data.seq{ismember(seq_data.gene,gene_list_measured(i))};
        length_gene(i,1) = length(seq);
    end
end

numerator = 0;
denominator = 0;
for i = 1:length(gene_list_measured)
    numerator = numerator + cost_gene(i)*length_gene(i)*abund_list_measured(i);
    denominator = denominator + length_gene(i)*abund_list_measured(i);
end

cost = numerator/denominator;


