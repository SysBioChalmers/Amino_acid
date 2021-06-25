load('cost_yeast_Yu_ClimR1.mat');
AA = cost_yeast.AA;
cost_gluc = cost_yeast.cost_gluc;
cost_limitedR1 = cost_yeast.cost_prot;
load('cost_yeast_Yu_ClimR2.mat');
cost_limitedR2 = cost_yeast.cost_prot;
load('cost_yeast_DiBartolomeo_GlucR1.mat');
cost_richR1 = cost_yeast.cost_prot;
load('cost_yeast_DiBartolomeo_GlucR2.mat');
cost_richR2 = cost_yeast.cost_prot;
load('cost_yeast_DiBartolomeo_GlucR3.mat');
cost_richR3 = cost_yeast.cost_prot;

cost = [cost_limitedR1 cost_limitedR2 cost_richR1 cost_richR2 cost_richR3 cost_gluc];
label = {'Glucose limited R1' 'Glucose limited R2' 'Glucose unlimited R1' 'Glucose unlimited R2' 'Glucose unlimited R3' 'Glucose cost'};

[RHO,~] = corr(cost,'type','Pearson');

figure();
h = heatmap(label,label,RHO,...
    'ColorMethod','count','CellLabelColor','k');
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';

figure('Name','PCA');
[~, score, ~, ~, explained, ~] = pca(cost','NumComponents',20);
h = scatter(score(:,1),score(:,2),30,'o','filled','LineWidth',1);
xlabel(['PC1 (',num2str(round(explained(1),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');
ylabel(['PC2 (',num2str(round(explained(2),1)),'%)'],'FontSize',7,'FontName','Helvetica','Color','k');


