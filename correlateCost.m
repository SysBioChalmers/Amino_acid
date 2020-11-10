% correlate cost

load('cost_yeast_DiBartolomeo_GlucR1.mat');
AA = cost_yeast.AA;
cost_gluc = cost_yeast.cost_gluc;
cost_DiBartolomeo1 = cost_yeast.cost_prot;
load('cost_yeast_DiBartolomeo_GlucR2.mat');
cost_DiBartolomeo2 = cost_yeast.cost_prot;
load('cost_yeast_DiBartolomeo_GlucR3.mat');
cost_DiBartolomeo3 = cost_yeast.cost_prot;
cost_DiBartolomeo = mean([cost_DiBartolomeo1 cost_DiBartolomeo2 cost_DiBartolomeo3],2);

load('cost_yeast_Yu_ClimR1.mat');
cost_Yu_Clim1 = cost_yeast.cost_prot;
load('cost_yeast_Yu_ClimR2.mat');
cost_Yu_Clim2 = cost_yeast.cost_prot;
cost_Yu_Clim = mean([cost_Yu_Clim1 cost_Yu_Clim2],2);

load('cost_yeast_Lahtvee_REF.mat');
cost_Lahtvee_REF = cost_yeast.cost_prot;

E_cost = [11.7;27.3;14.7;12.7;24.7
16.3;15.3;11.7;38.3;32.3
27.3;30.3;34.3;52;20.3
11.7;18.7;74.3;50;23.3];


% CplxScore = [4.76;56.34;33.72;32.72;57.16;
% 37.48;36.48;1;58.7;16.04;
% 16.04;30.14;64.68;44;31.8;
% 17.86;21.62;73;57;12.28]; % DOI: 10.1006/jtbi.1997.0443

%% cost_gluc vs e_cost
figure();
box on;
[RHO,~] = corr(cost_gluc,E_cost,'Type','Pearson');
text(cost_gluc,E_cost,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 3]);
ylim([0 80]);
text(0.1,75,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
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
title('batch','FontSize',6,'FontName','Helvetica');
xlabel('glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('protein cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 500 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

%% batch vs D0.2
figure();
box on;
[RHO,~] = corr(cost_DiBartolomeo,cost_Yu_Clim,'Type','Pearson');
text(cost_DiBartolomeo,cost_Yu_Clim,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 0.25]);
ylim([0 0.25]);
text(0.01,0.23,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost','FontSize',6,'FontName','Helvetica');
ylabel('glucose limited (D0.2)','FontSize',7,'FontName','Helvetica');
xlabel('glucose rich (batch)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 500 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);


%% batch vs D0.1
figure();
box on;
[RHO,~] = corr(cost_DiBartolomeo,cost_Lahtvee_REF,'Type','Pearson');
text(cost_DiBartolomeo,cost_Lahtvee_REF,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 0.9]);
ylim([0 0.9]);
text(0.25,0.83,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost','FontSize',6,'FontName','Helvetica');
ylabel('glucose limited (D0.1)','FontSize',7,'FontName','Helvetica');
xlabel('glucose rich (batch)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 600 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);




