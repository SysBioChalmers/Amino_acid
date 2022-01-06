load('YeastAA.mat');
AA = YeastAA.AAlist;

[~,txt,~] = xlsread(strcat('UniProt_Yeast.xlsx')); 
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);
data = zeros(length(AA),length(seq_data.gene));
rhodata = zeros(1,length(seq_data.gene));
pdata = zeros(1,length(seq_data.gene));
for i = 1:length(seq_data.gene)
    display([num2str(i) '/' num2str(length(seq_data.gene))]);
    seq = seq_data.seq{i};
    AAcount = zeros(length(AA),1);
    for j = 1:length(AA)
        AAcount(j,1) = length(strfind(seq,AA{j}));
    end
    data(:,i) = AAcount/sum(AAcount);
    [RHOtmp,PVALtmp] = corr(YeastAA.AAcontent_avg,AAcount/sum(AAcount),'Type','Pearson');
    rhodata(i) = RHOtmp;
    pdata(i) = PVALtmp;
end

proteinlist = seq_data.gene;
rholist = rhodata';
plist = pdata';


adjplist = mafdr(plist,'BHFDR', true);
