%% calculateCost4Protein
function cost_list = calculateCost4Protein(protein_list,AA,cost,org_name,tf_overlength)

[~,txt,~] = xlsread(strcat('UniProt_',org_name,'.xlsx')); 
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);

cost_list = zeros(length(protein_list),1);

for i = 1:length(protein_list)
    if ismember(protein_list(i),seq_data.gene)
        seq = seq_data.seq{ismember(seq_data.gene,protein_list(i))};
        AAcount = zeros(length(AA),1);
        for j = 1:length(AA)
            AAcount(j,1) = length(strfind(seq,AA{j}));
        end
        totcost = sum(AAcount.*cost);
        if tf_overlength == 1
            cost_list(i) = totcost/sum(AAcount);
        elseif tf_overlength == 0
            cost_list(i) = totcost;
        end
    end
end



