
expID = 'Yu_Clim';
% expID = 'DiBartolomeo_Gluc';

load(['cost_yeast_' expID '.mat']);

AA = cost_yeast.AA;
cost_1 = cost_yeast.cost_gluc;
cost_2 = cost_yeast.cost_prot;

[~,txt,~] = xlsread(strcat('UniProt_Yeast.xlsx')); 
seq_data.gene = txt(2:end,1);
seq_data.seq = txt(2:end,2);

data = zeros(length(AA),length(seq_data.gene));
for i = 1:length(seq_data.gene)
    display([num2str(i) '/' num2str(length(seq_data.gene))]);
    seq = seq_data.seq{i};
    AAcount = zeros(length(AA),1);
    for j = 1:length(AA)
        AAcount(j,1) = length(strfind(seq,AA{j}));
    end
    data(:,i) = AAcount/sum(AAcount);
    
end

% figure();
% dataheatmap = data;
% dataheatmap(dataheatmap>0.15) = 0.15;
% h = heatmap(1:length(seq_data.gene),AA,dataheatmap,'ColorMethod','count','CellLabelColor','none','GridVisible','off');
% h.Colormap = jet;
% set(h,'FontSize',6,'FontName','Helvetica');
% set(gcf,'position',[200 200 400 200]);
% set(gca,'position',[0.1 0.013 0.75 0.7]);


averageAA = calculateAAcomposition(seq_data.gene,AA,ones(length(seq_data.gene),1),'Yeast');

% figure();
% RHO = corr(averageAA,data,'Type','Pearson');
% [f,xi] = ksdensity(RHO);
% plot(xi,f);
% text(-0.5,3,['N = ' num2str(length(RHO))],'FontSize',7,'FontName','Helvetica');
% xlim([-0.6 1.1]);
% set(gca,'FontSize',6,'FontName','Helvetica');
% xlabel('Pearson r','FontSize',7,'FontName','Helvetica');
% ylabel('Density','FontSize',7,'FontName','Helvetica');
% set(gcf,'position',[600 400 100 60]);
% set(gca,'position',[0.2 0.35 0.7 0.55]);


figure();
subplot(2,1,1);
box on;
[RHO1,PVAL1] = corr(log10(averageAA),cost_1,'Type','Pearson');
text(log10(averageAA),cost_1,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',[49,130,189]/255);
xlim([-2.1 -0.9]);
ylim([0 3]);
text(-1.5,2.7,['Pearson r = ' num2str(round(RHO1,2))],'FontSize',6,'FontName','Helvetica');
text(-1.5,2.4,['p = ' num2str(round(PVAL1,10))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
ylabel('Glucose cost','FontSize',7,'FontName','Helvetica');

subplot(2,1,2);
box on;
[RHO4,PVAL4] = corr(log10(averageAA),cost_2,'Type','Pearson');
text(log10(averageAA),cost_2,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',[49,130,189]/255);
xlim([-2.1 -0.9]);
ylim([0 0.3]);
text(-1.5,0.27,['Pearson r = ' num2str(round(RHO4,2))],'FontSize',6,'FontName','Helvetica');
text(-1.5,0.24,['p = ' num2str(round(PVAL4,10))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
ylabel('Protein cost','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 140 220]);

figure();
hold on;
RHO1 = corr(cost_1,data,'Type','Pearson');
[f1,xi1] = ksdensity(RHO1);
plot(xi1,f1);
RHO4 = corr(cost_2,data,'Type','Pearson');
[f4,xi4] = ksdensity(RHO4);
plot(xi4,f4);
xlim([-1.1 0.6]);
legend('glucose cost','protein cost','Location','ne');
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Pearson r','FontSize',14,'FontName','Helvetica');
ylabel('Density','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[600 400 280 260]);
set(gca,'position',[0.15 0.15 0.8 0.8]);








