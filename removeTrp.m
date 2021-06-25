load('YeastAA.mat');

load('YeastCost/cost_kmax.mat');
AA = cost_kmax.AA;
cost_gluc = cost_kmax.cost_gluc;
cost_prot = cost_kmax.cost_prot;


E_cost = [14.5;20.5;18.5;15.5;26.5
10.5;9.5;14.5;29;38
37;36;36.5;61;14.5
14.5;21.5;75.5;59;29]; % ﻿Wagner (2005)

clr_gluc = [33,102,172]/255;
clr_prot = [178,24,43]/255;


idx = ~ismember(AA,'W');
AA = AA(idx);
cost_gluc = cost_gluc(idx);
cost_prot = cost_prot(idx);
E_cost = E_cost(idx);
YeastAA.AAcontent_avg = YeastAA.AAcontent_avg(idx);

figure();
subplot(1,3,2);
box on;
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
box on;
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
box on;
[RHO3,PVAL3] = corr(E_cost,log(YeastAA.AAcontent_avg*100),'Type','Pearson');
text(E_cost,log(YeastAA.AAcontent_avg*100),AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color','k');
xlim([0 90]);ylim([-0.5 2.5]);
title(['Pearson r = ' num2str(round(RHO3,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
ylabel('ln(% abundance)','FontSize',7,'FontName','Helvetica');
xlabel('Energy cost','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 300 80]);
