expList = {'Davidi_GLC_CHEM_mu=0.11_V'
'Davidi_GLC_CHEM_mu=0.12_S'
'Davidi_GLC_CHEM_mu=0.20_S'
'Davidi_GLC_CHEM_mu=0.21_P'
'Davidi_GLC_CHEM_mu=0.21_V'
'Davidi_GLC_CHEM_mu=0.22_P'
'Davidi_GLC_CHEM_mu=0.26_P'
'Davidi_GLC_CHEM_mu=0.31_P'
'Davidi_GLC_CHEM_mu=0.31_V'
'Davidi_GLC_CHEM_mu=0.35_S'
'Davidi_GLC_CHEM_mu=0.36_P'
'Davidi_GLC_CHEM_mu=0.40_V'
'Davidi_GLC_CHEM_mu=0.41_P'
'Davidi_GLC_CHEM_mu=0.46_P'
'Davidi_GLC_CHEM_mu=0.49_V'
'Davidi_GLC_CHEM_mu=0.50_S'
'Davidi_GLC_CHEM_mu=0.51_P'
'Davidi_GLC_BATCH_mu=0.58_S'};
cost_prot = zeros(20,0);
for i = 1:length(expList)
    load(['cost_ecoli_' expList{i} '.mat']);
    cost_prot = [cost_prot cost_ecoli.cost_prot];
end
AA = cost_ecoli.AA;
cost_gluc = cost_ecoli.cost_gluc;

E_cost = [11.7;27.3;14.7;12.7;24.7
16.3;15.3;11.7;38.3;32.3
27.3;30.3;34.3;52;20.3
11.7;18.7;74.3;50;23.3];

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

figure();
RHO = corr(cost_prot);
heatmap(RHO)

figure();
box on;
a = 1;
b = 18;
[RHO,~] = corr(cost_prot(:,a),cost_prot(:,b),'Type','Pearson');
text(cost_prot(:,a),cost_prot(:,b),AA,'VerticalAlignment','middle','HorizontalAlignment','center','FontSize',6,'FontName','Helvetica');
xlim([0 0.95]);
ylim([0 0.95]);
text(0.1,0.87,['Pearson r = ' num2str(round(RHO,2))],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
title('protein cost','FontSize',6,'FontName','Helvetica');
xlabel(['condition' num2str(a)],'FontSize',7,'FontName','Helvetica');
ylabel(['condition' num2str(b)],'FontSize',7,'FontName','Helvetica');
set(gcf,'position',[800 400 100 100]);
set(gca,'position',[0.25 0.25 0.68 0.68]);

