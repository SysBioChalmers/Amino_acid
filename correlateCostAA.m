load('YeastAA.mat');
n = size(YeastAA.AAcontent,2);
rhodata = zeros(n,n);
pdata = zeros(n,n);
for i = 1:n
    idata = YeastAA.AAcontent(:,i);
    for j = 1:n
        jdata = YeastAA.AAcontent(:,j);
        [RHOtmp,PVALtmp] = corr(idata,jdata,'Type','Pearson');
        rhodata(i,j) = RHOtmp;
        pdata(i,j) = PVALtmp;
    end
end

figure();
xlbl = num2cell(1:1:n);
ylbl = num2cell(1:1:n);
minclr = [255,237,160]/255;
maxclr = [240,59,32]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(xlbl,ylbl,rhodata,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title('Correlations between amino acid relative abundances');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[820 320 360 360]);
set(gca,'position',[0.25 0.2 0.5 0.5]);

figure();
% box on;
[RHO,~] = corr(YeastAA.AAcontent(:,end)*100,YeastAA.AAcontent(:,1)*100,'Type','Pearson');
text(YeastAA.AAcontent(:,end)*100,YeastAA.AAcontent(:,1)*100,YeastAA.AAlist,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 10]);
ylim([0 10]);
% text(0.1,73,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('% abundance','FontSize',6,'FontName','Helvetica');
xlabel('Condition 30','FontSize',7,'FontName','Helvetica');
ylabel('Condition 1','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 80 80]);
set(gca,'position',[0.27 0.27 0.62 0.62]);



cd YeastCost/;
load('cost_kmax.mat');
AA = cost_kmax.AA;
cost_gluc = cost_kmax.cost_gluc;
cost_prot = cost_kmax.cost_prot;
cd ../;

% E_cost = [11.7;27.3;14.7;12.7;24.7
% 16.3;15.3;11.7;38.3;32.3
% 27.3;30.3;34.3;52;20.3
% 11.7;18.7;74.3;50;23.3];

E_cost = [14.5;20.5;18.5;15.5;26.5
10.5;9.5;14.5;29;38
37;36;36.5;61;14.5
14.5;21.5;75.5;59;29]; % ﻿Wagner (2005)

clr_gluc = [33,102,172]/255;
clr_prot = [178,24,43]/255;

figure();
subplot(1,3,2);
% box on;
[RHO1,PVAL1] = corr(cost_gluc,log(YeastAA.AAcontent_avg*100),'Type','Pearson');
text(cost_gluc,log(YeastAA.AAcontent_avg*100),AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_gluc);
xlim([0 2.5]);ylim([-0.5 2.5]);
% text(1.5,2,['Pearson r = ' num2str(round(RHO1,2))],'FontSize',6,'FontName','Helvetica');
% text(2,2,['p = ' num2str(round(PVAL1,3))],'FontSize',6,'FontName','Helvetica');
title(['Pearson r = ' num2str(round(RHO1,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('ln(% abundance)','FontSize',7,'FontName','Helvetica');
xlabel('Glucose cost','FontSize',7,'FontName','Helvetica');

subplot(1,3,3);
% box on;
[RHO2,PVAL2] = corr(cost_prot,log(YeastAA.AAcontent_avg*100),'Type','Pearson');
text(cost_prot,log(YeastAA.AAcontent_avg*100),AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_prot);
xlim([0 0.25]);ylim([-0.5 2.5]);
% text(0.15,2,['Pearson r = ' num2str(round(RHO2,2))],'FontSize',6,'FontName','Helvetica');
% text(0.2,2,['p = ' num2str(round(PVAL2,8))],'FontSize',6,'FontName','Helvetica');
title(['Pearson r = ' num2str(round(RHO2,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('ln(% abundance)','FontSize',7,'FontName','Helvetica');
xlabel('Protein cost (kmax)','FontSize',7,'FontName','Helvetica');

subplot(1,3,1);
% box on;
[RHO3,PVAL3] = corr(E_cost,log(YeastAA.AAcontent_avg*100),'Type','Pearson');
text(E_cost,log(YeastAA.AAcontent_avg*100),AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color','k');
xlim([0 90]);ylim([-0.5 2.5]);
title(['Pearson r = ' num2str(round(RHO3,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('ln(% abundance)','FontSize',7,'FontName','Helvetica');
xlabel('Energy cost','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 300 80]);



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
    [RHO1,PVAL1] = corr(log(AA_tmp*100),cost_kapp.cost_gluc,'Type','Pearson');
    [RHO2,PVAL2] = corr(log(AA_tmp*100),cost_kapp.cost_prot,'Type','Pearson');
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





