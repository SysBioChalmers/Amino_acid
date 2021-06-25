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