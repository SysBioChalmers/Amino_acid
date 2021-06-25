% import molecular weight
[num,txt,~] = xlsread(strcat('UniProt_Yeast.xlsx'));
yeast_protein = txt(2:end,1);
yeast_mw = num;
clear txt num;

YeastAAall.condition = cell(1,0);
YeastAAall.AAlist = {'A';'R';'N';'D';'C';
    'Q';'E';'G';'H';'I';
    'L';'K';'M';'F';'P';
    'S';'T';'W';'Y';'V'};
YeastAAall.AAcontent = zeros(length(YeastAAall.AAlist),0);

[num,txt,~] = xlsread('../Yeast_kapp/kappEstimation/Data/ProteomicsFlux.xlsx','Lahtvee2017');
condlist = txt(1,2:end);
protlist = txt(2:end,1);
for i = 1:length(condlist)
    display([num2str(i),'/',num2str(length(condlist))]);
    YeastAAall.condition = [YeastAAall.condition condlist(i)];
    AAcomposition = calculateAAcomposition(protlist,YeastAAall.AAlist,num(:,i),'Yeast');
    YeastAAall.AAcontent = [YeastAAall.AAcontent AAcomposition];
end
clear txt num;

[num,txt,~] = xlsread('../Yeast_kapp/kappEstimation/Data/ProteomicsFlux.xlsx','Yu2020');
condlist = txt(1,2:end);
protlist = txt(2:end,1);
for i = 1:length(condlist)
    display([num2str(i),'/',num2str(length(condlist))]);
    YeastAAall.condition = [YeastAAall.condition condlist(i)];
    AAcomposition = calculateAAcomposition(protlist,YeastAAall.AAlist,num(:,i),'Yeast');
    YeastAAall.AAcontent = [YeastAAall.AAcontent AAcomposition];
end
clear txt num;

[num,txt,~] = xlsread('../Yeast_kapp/kappEstimation/Data/ProteomicsFlux.xlsx','DiBartolomeo2020');
condlist = txt(1,2:end);
protlist = txt(2:end,1);
for i = 1:length(condlist)
    display([num2str(i),'/',num2str(length(condlist))]);
    YeastAAall.condition = [YeastAAall.condition condlist(i)];
    
    abund_tmp = num(:,i);
    [~,b] = ismember(protlist,yeast_protein);
    protlist = protlist(b~=0);
    abund_tmp = abund_tmp(b~=0);
    b = b(b~=0);
    mwlist = yeast_mw(b);
    
    AAcomposition = calculateAAcomposition(protlist,YeastAAall.AAlist,abund_tmp./mwlist,'Yeast');
    YeastAAall.AAcontent = [YeastAAall.AAcontent AAcomposition];
end
clear txt num;

[num,txt,~] = xlsread('../Yeast_kapp/kappEstimation/Data/ProteomicsFlux.xlsx','Yu2021');
condlist = txt(1,2:end);
protlist = txt(2:end,1);
for i = 1:length(condlist)
    display([num2str(i),'/',num2str(length(condlist))]);
    YeastAAall.condition = [YeastAAall.condition condlist(i)];
    AAcomposition = calculateAAcomposition(protlist,YeastAAall.AAlist,num(:,i),'Yeast');
    YeastAAall.AAcontent = [YeastAAall.AAcontent AAcomposition];
end
clear txt num;


conditions = strrep(YeastAAall.condition,'R1','');
conditions = strrep(conditions,'R2','');
conditions = strrep(conditions,'R3','');
conditions = unique(conditions);
YeastAA = struct();
YeastAA.AAlist = YeastAAall.AAlist;
YeastAA.condition = conditions;
YeastAA.AAcontent = zeros(length(YeastAA.AAlist),length(conditions));
for i = 1:length(YeastAA.condition)
    idx = contains(YeastAAall.condition,YeastAA.condition(i));
    YeastAA.AAcontent(:,i) = mean(YeastAAall.AAcontent(:,idx),2);
end

YeastAA.AAcontent_avg = mean(YeastAA.AAcontent,2);
save('YeastAA.mat','YeastAA');





