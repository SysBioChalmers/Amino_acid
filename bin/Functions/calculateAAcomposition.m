%% calculateAAcomposition
function composition_list = calculateAAcomposition(gene_list,AA_list,abund_list,org_name)

[~,txt,~] = xlsread(strcat('UniProt_',org_name,'.xlsx')); 
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);

abs = zeros(length(AA_list),1);
for i = 1:length(gene_list)
    if ismember(gene_list(i),seq_data.gene)
        seq = seq_data.seq{ismember(seq_data.gene,gene_list(i))};
        AAcount = zeros(length(AA_list),1);
        for j = 1:length(AA_list)
            AAcount(j,1) = length(strfind(seq,AA_list{j}));
        end
        abs = abs + AAcount*abund_list(i);
    end
end

% composition_list = abs;

composition_list = abs/sum(abs);

