% correlate cost

load('cost_yeast_DiBartolomeo_Gluc.mat');
AA = cost_yeast.AA;
cost_gluc = cost_yeast.cost_gluc;
cost_DiBartolomeo = cost_yeast.cost_prot;

load('cost_yeast_Yu_Clim.mat');
cost_Yu_Clim = cost_yeast.cost_prot;

E_cost = [11.7;27.3;14.7;12.7;24.7
16.3;15.3;11.7;38.3;32.3
27.3;30.3;34.3;52;20.3
11.7;18.7;74.3;50;23.3];

% CplxScore = [4.76;56.34;33.72;32.72;57.16;
% 37.48;36.48;1;58.7;16.04;
% 16.04;30.14;64.68;44;31.8;
% 17.86;21.62;73;57;12.28]; % DOI: 10.1006/jtbi.1997.0443

% N_atoms = [1;4;2;1;1;
% 2;1;1;3;1;
% 1;2;1;1;1;
% 1;1;2;1;1]; % DOI: 10.1126/SCIENCE.AAZ9642

%% cost_gluc vs e_cost
figure();
box on;
[RHO,~] = corr(cost_gluc,E_cost,'Type','Pearson');
text(cost_gluc,E_cost,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 3]);
ylim([0 80]);
text(0.1,73,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('energy cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% batch
figure();
box on;
[RHO,~] = corr(cost_gluc,cost_DiBartolomeo,'Type','Pearson');
text(cost_gluc,cost_DiBartolomeo,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
ylim([0 0.25]);
xlim([0 3]);
text(0.1,0.23,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
% title('batch','FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('protein cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[300 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% batch vs D0.2
figure();
box on;
[RHO,~] = corr(cost_DiBartolomeo,cost_Yu_Clim,'Type','Pearson');
text(cost_DiBartolomeo,cost_Yu_Clim,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 0.35]);
ylim([0 0.35]);
text(0.01,0.32,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost','FontSize',6,'FontName','Helvetica');
ylabel('glucose limited (chemostat)','FontSize',7,'FontName','Helvetica');
xlabel('glucose rich (batch)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 500 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);



%% ferment

load('ferment_cost_yeast_DiBartolomeo_Gluc.mat');
cost_gluc_f = cost_yeast.cost_gluc;
cost_DiBartolomeo_f = cost_yeast.cost_prot;

figure();
box on;
[RHO,~] = corr(cost_DiBartolomeo,cost_DiBartolomeo_f,'Type','Pearson');
text(cost_DiBartolomeo,cost_DiBartolomeo_f,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 0.49]);
ylim([0 0.49]);
text(0.15,0.45,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost (batch kapp)','FontSize',6,'FontName','Helvetica');
xlabel('respiration (default)','FontSize',7,'FontName','Helvetica');
ylabel('fermentation','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[100 300 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% ethanol as CS
load('ethanol_cost_yeast_DiBartolomeo_Etoh.mat');
cost_etoh = cost_yeast.cost_gluc;
cost_DiBartolomeo_e = cost_yeast.cost_prot;

figure();
box on;
[RHO,~] = corr(cost_gluc,cost_etoh,'Type','Pearson');
text(cost_gluc,cost_etoh,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 8]);
ylim([0 8]);
text(2.5,0.7,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('ethanol cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 500 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

figure();
box on;
[RHO,~] = corr(cost_DiBartolomeo,cost_DiBartolomeo_e,'Type','Pearson');
text(cost_DiBartolomeo,cost_DiBartolomeo_e,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 1.15]);
ylim([0 1.15]);
text(0.35,0.1,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost','FontSize',6,'FontName','Helvetica');
xlabel('glucose as CS','FontSize',7,'FontName','Helvetica');
ylabel('ethanol as CS','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 600 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);








