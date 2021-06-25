cd YeastCost/;
load('cost_kmax.mat');
AA = cost_kmax.AA;
cost_gluc = cost_kmax.cost_gluc;
cost_prot_kmax = cost_kmax.cost_prot;

file = dir('cost_kapp_*.mat');
cost_prot_kapp = zeros(length(AA),length(file));
cond_name = cell(1,length(file));
for i = 1:length(file)
    display([num2str(i) '/' num2str(length(file))]);
    filename = file(i).name;
    load(filename);
    cost_prot_kapp(:,i) = cost_kapp.cost_prot;
    condname = strrep(filename,'cost_kapp_','');
    condname = strrep(condname,'.mat','');
    cond_name(i) = {condname};
end

cd ../;

cost = [cost_gluc cost_prot_kmax cost_prot_kapp];

n = size(cost_prot_kapp,2)+2;
rhodata = zeros(n,n);
pdata = zeros(n,n);
for i = 1:n
    idata = cost(:,i);
    for j = 1:n
        jdata = cost(:,j);
        [RHOtmp,PVALtmp] = corr(idata,jdata,'Type','Pearson');
        rhodata(i,j) = RHOtmp;
        pdata(i,j) = PVALtmp;
    end
end

figure();
xlbl = ['Glucose cost' 'Protein cost (kmax)' cond_name];
ylbl = xlbl;
minclr = [255,237,160]/255;
maxclr = [240,59,32]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(xlbl,ylbl,rhodata,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title('Correlations between amino acid costs');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[820 320 360 360]);
set(gca,'position',[0.25 0.2 0.5 0.5]);

%% cost_gluc vs e_cost
% E_cost = [11.7;27.3;14.7;12.7;24.7
% 16.3;15.3;11.7;38.3;32.3
% 27.3;30.3;34.3;52;20.3
% 11.7;18.7;74.3;50;23.3];

E_cost = [14.5;20.5;18.5;15.5;26.5
10.5;9.5;14.5;29;38
37;36;36.5;61;14.5
14.5;21.5;75.5;59;29]; % ﻿Wagner (2005)

figure();
box on;
[RHO,~] = corr(cost_gluc,E_cost,'Type','Pearson');
text(cost_gluc,E_cost,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 3]);
ylim([0 80]);
% text(0.1,73,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title(['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
xlabel('Glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('Energy cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 100 100]);
set(gca,'position',[0.27 0.2 0.68 0.68]);

%% cost_gluc vs cost_prot (kmax)
figure();
box on;
[RHO,~] = corr(cost_gluc,cost_prot_kmax,'Type','Pearson');
text(cost_gluc,cost_prot_kmax,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
ylim([-0.05 0.25]);
xlim([0 3]);
% text(0.1,0.23,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title(['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
xlabel('Glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('Protein cost (kmax)','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[300 400 100 100]);
set(gca,'position',[0.27 0.2 0.68 0.68]);




