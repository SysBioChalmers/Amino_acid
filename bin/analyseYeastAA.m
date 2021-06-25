
% expID = 'Yu_Clim';
expID = 'DiBartolomeo_Gluc';

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

figure();
dataheatmap = data;
dataheatmap(dataheatmap>0.15) = 0.15;
h = heatmap(1:length(seq_data.gene),AA,dataheatmap,'ColorMethod','count','CellLabelColor','none','GridVisible','off');
h.Colormap = jet;
set(h,'FontSize',6,'FontName','Helvetica');
set(gcf,'position',[200 200 400 200]);
set(gca,'position',[0.1 0.013 0.75 0.7]);


averageAA = calculateAAcomposition(seq_data.gene,AA,ones(length(seq_data.gene),1),'Yeast');

figure();
RHO = corr(averageAA,data,'Type','Pearson');
[f,xi] = ksdensity(RHO);
plot(xi,f);
text(-0.5,3,['N = ' num2str(length(RHO))],'FontSize',7,'FontName','Helvetica');
xlim([-0.6 1.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Pearson r','FontSize',7,'FontName','Helvetica');
ylabel('Density','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 100 60]);
set(gca,'position',[0.2 0.35 0.7 0.55]);


clr_gluc = [5,113,176]/255;
clr_prot = [202,0,32]/255;


figure();
subplot(2,1,1);
box on;
[RHO1,PVAL1] = corr(log10(averageAA),cost_1,'Type','Pearson');
text(log10(averageAA),cost_1,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_gluc);
xlim([-2.1 -0.9]);
ylim([0 2.45]);
text(-2.05,0.5,['Pearson r = ' num2str(round(RHO1,2))],'FontSize',6,'FontName','Helvetica');
text(-2.05,0.2,['p = ' num2str(round(PVAL1,3))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
ylabel('Glucose cost','FontSize',7,'FontName','Helvetica');

subplot(2,1,2);
box on;
[RHO2,PVAL2] = corr(log10(averageAA),cost_2,'Type','Pearson');
text(log10(averageAA),cost_2,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_prot);
xlim([-2.1 -0.9]);
ylim([0 0.245]);
text(-2.05,0.05,['Pearson r = ' num2str(round(RHO2,2))],'FontSize',6,'FontName','Helvetica');
text(-2.05,0.02,['p = ' num2str(round(PVAL2,6))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
ylabel('Protein cost','FontSize',7,'FontName','Helvetica');

set(gcf,'position',[200 200 140 220]);

figure();
hold on;
RHO1 = corr(cost_1,data,'Type','Pearson');
[f1,xi1] = ksdensity(RHO1);
plot(xi1,f1,'Color',clr_gluc);
RHO2 = corr(cost_2,data,'Type','Pearson');
[f2,xi2] = ksdensity(RHO2);
plot(xi2,f2,'Color',clr_prot);
xlim([-1.1 0.6]);
legend('glucose cost','protein cost','Location','ne');
set(gca,'FontSize',12,'FontName','Helvetica');
xlabel('Pearson r','FontSize',14,'FontName','Helvetica');
ylabel('Density','FontSize',14,'FontName','Helvetica');
set(gcf,'position',[600 400 280 260]);
set(gca,'position',[0.15 0.15 0.8 0.8]);








