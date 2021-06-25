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