% load('GEM-ecoli.mat');
% model_split = splitRevRxns(model);
% save('GEM-ecoli-split.mat','model_split');

load('GEM-ecoli-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'escherichia coli');

enzymedata = updateEcoliEnzyme(enzymedata,model_split);

enzymedata_org = enzymedata;

step = 0.001;
mu = 1;
tol = mu*0.00001;
substrateExID = 'EX_glc__D_e'; % glucose

load('EcoliSA.mat');

cond = EcoliSA.condition(contains(EcoliSA.condition,'GLC','IgnoreCase',true));

prot_cost = zeros(20,length(cond));

for k = 1:length(cond)
    display([num2str(k),'/',num2str(length(cond))]);
    enzymedata = enzymedata_org;
    enzymedata = updateCondSpecEcoliEnzyme(enzymedata,cond(k));
    
    model = convertModel(model_split,enzymedata);
    
    model = blockEcoliRxns(model);
    
    
    cost_ecoli = calculateEcoliCost(model,substrateExID,mu,step,tol);
    
    
    AA = cost_ecoli.AA;
    cost_1 = cost_ecoli.cost_gluc;
    cost_2 = cost_ecoli.cost_prot;
    
    prot_cost(:,k) = cost_2;
    
%     [~,txt,~] = xlsread(strcat('UniProt_Ecoli.xlsx'));
%     seq_data.gene = txt(2:end,1);
%     seq_data.seq = txt(2:end,2);
%     
%     averageAA = calculateAAcomposition(seq_data.gene,AA,ones(length(seq_data.gene),1),'Ecoli');
%     
%     clr_gluc = [5,113,176]/255;
%     clr_prot = [202,0,32]/255;
%     
%     figure();
%     subplot(2,1,1);
%     box on;
%     [RHO1,PVAL1] = corr(log10(averageAA),cost_1,'Type','Pearson');
%     text(log10(averageAA),cost_1,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_gluc);
%     xlim([-2.1 -0.9]);
%     ylim([0 2.45]);
%     text(-2.05,0.5,['Pearson r = ' num2str(round(RHO1,2))],'FontSize',6,'FontName','Helvetica');
%     text(-2.05,0.2,['p = ' num2str(round(PVAL1,3))],'FontSize',6,'FontName','Helvetica');
%     set(gca,'FontSize',6,'FontName','Helvetica');
%     xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
%     ylabel('Glucose cost','FontSize',7,'FontName','Helvetica');
%     
%     subplot(2,1,2);
%     box on;
%     [RHO2,PVAL2] = corr(log10(averageAA),cost_2,'Type','Pearson');
%     text(log10(averageAA),cost_2,AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica','Color',clr_prot);
%     xlim([-2.1 -0.9]);
%     ylim([0 0.245]);
%     text(-2.05,0.05,['Pearson r = ' num2str(round(RHO2,2))],'FontSize',6,'FontName','Helvetica');
%     text(-2.05,0.02,['p = ' num2str(round(PVAL2,6))],'FontSize',6,'FontName','Helvetica');
%     set(gca,'FontSize',6,'FontName','Helvetica');
%     xlabel('log10(Average AA)','FontSize',7,'FontName','Helvetica');
%     ylabel('Protein cost','FontSize',7,'FontName','Helvetica');
%     
%     set(gcf,'position',[200 200 140 220]);
    
end

cost = [cost_1 prot_cost];

n = length(cond)+1;
rhodata = zeros(n,n);
pdata = zeros(n,n);
for i = 1:n
    idata = cost(:,i);
    for j = 1:n
        jdata = cost(:,j);
        [RHOtmp,PVALtmp] = corr(idata,jdata,'Type','Spearman');
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
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[820 320 360 360]);
set(gca,'position',[0.25 0.2 0.5 0.5]);


