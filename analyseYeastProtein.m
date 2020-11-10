load('cost_min_gluc_yeast.mat');
load('cost_min_prot_yeast.mat');

AA = cost_min_gluc_yeast.AA;
cost_1 = cost_min_gluc_yeast.cost_gluc;
cost_4 = cost_min_prot_yeast.cost_prot;

[~,txt,~] = xlsread(strcat('UniProt_Yeast.xlsx')); 
gene_list = txt(2:end,1);

averageAA = calculateAAcomposition(gene_list,AA,ones(length(gene_list),1),'Yeast');

[~, cost_average1] = costPerProtein(gene_list,AA,cost_1,'Yeast');
[~, cost_average4] = costPerProtein(gene_list,AA,cost_4,'Yeast');

averagecost_1 = sum(averageAA.*cost_1);
averagecost_4 = sum(averageAA.*cost_4);

% dv = 0.1;
% q = (cost_average1/averagecost_1)./(cost_average4/averagecost_4);
% idx = q < 1+dv & q > 1-dv;
% gene_list = gene_list(~idx);
% cost_average1 = cost_average1(~idx);
% cost_average4 = cost_average4(~idx);

figure();
hold on;
box on;
line([0 3],[0 3],'Color','k');
scatter(cost_average1/averagecost_1,cost_average4/averagecost_4,5,'filled','MarkerEdgeColor',[.5 .5 .5],...
              'MarkerFaceColor',[.7 .7 .7],'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3,'LineWidth',0.5);
xlim([0.4 2.1]);
ylim([0.4 2.1]);
set(gca,'FontSize',6,'FontName','Helvetica');
text(1.3,0.6,['N = ' num2str(length(gene_list))],'FontSize',7,'FontName','Helvetica');
xlabel('Normalized glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('Normalized protein cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);





% [pax_abund,txt,~] = xlsread('paxdb_yeast.xlsx');
% pax_genes = txt(2:end,1);
% clear txt;
% [M_idx,~] = ismember(gene_list,pax_genes);
% pax_gene_list = gene_list(M_idx);
% pax_abund = pax_abund(M_idx,:);
% [~, pax_cost_average1] = costPerProtein(pax_gene_list,AA,cost_1,'Yeast');
% [~, pax_cost_average4] = costPerProtein(pax_gene_list,AA,cost_4,'Yeast');
% figure();
% subplot(2,1,1);
% scatter(log10(pax_abund),pax_cost_average1,5,'filled','MarkerFaceAlpha',0.3);
% xlim([-4 6]);
% % ylim([0.4 2.1]);
% set(gca,'FontSize',6,'FontName','Helvetica');
% text(-3.7,1.35,['N = ' num2str(length(pax_gene_list))],'FontSize',7,'FontName','Helvetica');
% xlabel('log10(abundance)','FontSize',7,'FontName','Helvetica');
% ylabel('Glucose cost','FontSize',7,'FontName','Helvetica');
% subplot(2,1,2);
% scatter(log10(pax_abund),pax_cost_average4,5,'filled','MarkerFaceAlpha',0.3);
% xlim([-4 6]);
% % ylim([0.4 2.1]);
% set(gca,'FontSize',6,'FontName','Helvetica');
% text(-3.7,1.35,['N = ' num2str(length(pax_gene_list))],'FontSize',7,'FontName','Helvetica');
% xlabel('log10(abundance)','FontSize',7,'FontName','Helvetica');
% ylabel('Protein cost','FontSize',7,'FontName','Helvetica');
% set(gcf,'position',[600 200 100 200]);



% use fc 
[data_abund,txt] = xlsread('proteomics_data_yeast.xlsx','Sheet_protein'); % unpublished
gene_list = txt(2:end,1);
label_list = txt(1,2:end);
clear txt;
fc = data_abund(:,9)./data_abund(:,3);
idxnan = isnan(fc);
gene_list = gene_list(~idxnan);
fc = fc(~idxnan);
log2fc = log2(fc);


[~, cost_average1] = costPerProtein(gene_list,AA,cost_1,'Yeast');
[~, cost_average4] = costPerProtein(gene_list,AA,cost_4,'Yeast');


% scatter(log2fc,cost_average1./cost_average4,10,'filled');

figure();
hold on;
box on;
line([0 3],[0 3],'Color','k');
scatter(cost_average1/averagecost_1,cost_average4/averagecost_4,10,log2fc,'filled');
xlim([0.5 1.5]);
ylim([0.5 1.5]);
set(gca,'FontSize',6,'FontName','Helvetica');
text(1,0.6,['N = ' num2str(length(gene_list))],'FontSize',7,'FontName','Helvetica');
colorbar;
xlabel('Normalized glucose cost','FontSize',7,'FontName','Helvetica');
ylabel('Normalized protein cost','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 150 100]);
set(gca,'position',[0.25 0.25 0.48 0.68]);

% x = (cost_average1/averagecost_1)./(cost_average4/averagecost_4);
% y = log2fc;
% scatter(x,y);




