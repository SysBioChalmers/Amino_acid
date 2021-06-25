clr_gluc = [33,102,172]/255;
clr_prot = [178,24,43]/255;

load('YeastAA.mat');

cd YeastCost/;

file = dir('cost_kapp_*.mat');
rho = zeros(0,2);
pval = zeros(0,2);
condition = cell(0,1);
for i = 1:length(file)
    display([num2str(i) '/' num2str(length(file))]);
    filename = file(i).name;
    load(filename);
    condID = strrep(filename,'cost_kapp_','');
    condID = strrep(condID,'.mat','');
    AA_tmp = YeastAA.AAcontent(:,ismember(YeastAA.condition,condID));
    [RHO1,PVAL1] = corr(log(AA_tmp),cost_kapp.cost_gluc,'Type','Pearson');
    [RHO2,PVAL2] = corr(log(AA_tmp),cost_kapp.cost_prot,'Type','Pearson');
    condition = [condition;{condID}];
    rho = [rho;[RHO1 RHO2]];
    pval = [pval;[PVAL1 PVAL2]];
end
cd ../;

figure();
b = bar(1:length(condition),rho,0.7);
b(1).LineWidth = 0.01;
b(1).FaceColor = clr_gluc;
b(2).LineWidth = 0.01;
b(2).FaceColor = clr_prot;
set(gca,'XTick',1:1:length(condition));
set(gca,'XTickLabel',num2cell(1:1:length(condition)));
legend({'Glucose cost' 'Protein cost'},'Location','northwest','Orientation','horizontal','FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('Pearson r','FontSize',7,'FontName','Helvetica','Color','k');
title('Correlations between amino acid costs and relative abundances','FontSize',7,'FontName','Helvetica','Color','k');
set(gcf,'position',[200 100 300 125]);
set(gca,'position',[0.1 0.2 0.75 0.5]);
% box off;