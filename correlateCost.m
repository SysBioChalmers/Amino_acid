% correlate cost

load('cost_min_gluc_yeast.mat');
load('cost_min_prot_yeast.mat');

AA = cost_min_gluc_yeast.AA;

E_cost = [11.7;27.3;14.7;12.7;24.7
16.3;15.3;11.7;38.3;32.3
27.3;30.3;34.3;52;20.3
11.7;18.7;74.3;50;23.3];

cost_1 = cost_min_gluc_yeast.cost_gluc;
cost_2 = cost_min_gluc_yeast.cost_prot;
cost_3 = cost_min_prot_yeast.cost_gluc;
cost_4 = cost_min_prot_yeast.cost_prot;

% CplxScore = [4.76;56.34;33.72;32.72;57.16;
% 37.48;36.48;1;58.7;16.04;
% 16.04;30.14;64.68;44;31.8;
% 17.86;21.62;73;57;12.28]; % DOI: 10.1006/jtbi.1997.0443

%% gluc_gluc vs gluc_prot
figure();
box on;
scatter(cost_1,cost_2,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_1,cost_2,'Type','Pearson');
% xlim([-1 4]);
% ylim([-1 4]);
text(0.1,0.14,['Pearson r = ' num2str(round(RHO,3))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('cost 1','FontSize',7,'FontName','Helvetica');
ylabel('cost 2','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% prot_gluc vs prot_prot
figure();
box on;
scatter(cost_3,cost_4,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_3,cost_4,'Type','Pearson');
% xlim([-1 4]);
% ylim([-1 4]);
text(0.33,0.14,['Pearson r = ' num2str(round(RHO,3))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('cost 3','FontSize',7,'FontName','Helvetica');
ylabel('cost 4','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% gluc_gluc vs gluc_gluc
figure();
box on;
scatter(cost_1,cost_3,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_1,cost_3,'Type','Pearson');
% xlim([-1 4]);
% ylim([-1 4]);
text(0.9,0.8,['Pearson r = ' num2str(round(RHO,3))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('cost 1','FontSize',7,'FontName','Helvetica');
ylabel('cost 3','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[400 200 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% prot_prot vs prot_prot
figure();
box on;
scatter(cost_2,cost_4,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_2,cost_4,'Type','Pearson');
% xlim([-1 4]);
% ylim([-1 4]);
text(0.045,0.012,['Pearson r = ' num2str(round(RHO,3))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('cost 2','FontSize',7,'FontName','Helvetica');
ylabel('cost 4','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[400 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% gluc_gluc vs prot_prot
figure();
box on;
% scatter(cost_1,cost_4,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_1,cost_4,'Type','Pearson');
text(cost_1,cost_4,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 3]);
ylim([0 0.15]);
text(0.1,0.14,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('protein cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 600 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% gluc_gluc vs e_cost
figure();
box on;
% scatter(cost_1,E_cost,10,'o','filled','LineWidth',1,'MarkerEdgeColor',[49,130,189]/255,'MarkerFaceColor',[49,130,189]/255,'MarkerFaceAlpha',0.5);
[RHO,~] = corr(cost_1,E_cost,'Type','Pearson');
text(cost_1,E_cost,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 3]);
ylim([0 80]);
text(0.1,75,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('energy cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);



