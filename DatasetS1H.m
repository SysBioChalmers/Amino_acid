
load('YeastCost/cost_kmax.mat');

load('GEM-yeast-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);
model = convertModel(model_split,enzymedata);
clear model_split;

pcf = zeros(length(cost_kmax.AA),6);
for i = 1:length(cost_kmax.AA)
    disp([num2str(i),'/',num2str(length(cost_kmax.AA))]);
    fluxAA = cost_kmax.flux_AAs(:,i);
    
    flux_tmp = fluxAA - cost_kmax.flux_ref;
    active_rxns = model.rxns(flux_tmp ~= 0);
    active_fluxes = flux_tmp(flux_tmp ~= 0);
    protein_cost = zeros(length(active_rxns),1);
    protein_cost_cs = zeros(length(active_rxns),1);
    
    for j = 1:length(active_rxns)
        if ismember(active_rxns(j),enzymedata.rxn)
            protein_cost(j) = enzymedata.minPC(ismember(enzymedata.rxn,active_rxns(j)));
            protein_cost_cs(j) = enzymedata.kcat_conf(ismember(enzymedata.rxn,active_rxns(j)));
        end
    end
    totAbsPC = sum(abs(active_fluxes.*protein_cost));
    
    AbsPC5 = sum(abs(protein_cost(protein_cost_cs==5).*active_fluxes(protein_cost_cs==5)));
    AbsPC4 = sum(abs(protein_cost(protein_cost_cs==4).*active_fluxes(protein_cost_cs==4)));
    AbsPC3 = sum(abs(protein_cost(protein_cost_cs==3).*active_fluxes(protein_cost_cs==3)));
    AbsPC2 = sum(abs(protein_cost(protein_cost_cs==2).*active_fluxes(protein_cost_cs==2)));
    AbsPC1 = sum(abs(protein_cost(protein_cost_cs==1).*active_fluxes(protein_cost_cs==1)));
    AbsPC0 = sum(abs(protein_cost(protein_cost_cs==0).*active_fluxes(protein_cost_cs==0)));
    
    pcf(i,:) = [AbsPC5 AbsPC4 AbsPC3 AbsPC2 AbsPC1 AbsPC0] / totAbsPC;
    
end

data = [pcf(:,1) sum(pcf(:,2:3),2) sum(pcf(:,4:6),2)];

%% 

figure();
hold on;

databar = data*100;

b = bar(databar,'stacked');

b(3).FaceColor = [253,219,199]/255;
b(2).FaceColor = [214,96,77]/255;
b(1).FaceColor = [178,24,43]/255;

xticks(1:1:20);
xticklabels(cost_kmax.AA);
xlim([0 21]);
ylim([0 100]);
ylabel('Contribution to total protein cost (%)','FontSize',7,'FontName','Helvetica');
set(gca,'FontSize',7,'FontName','Helvetica');
legend({'in vivo kmax' 'yeast kcat' 'non-yeast kcat'},'FontSize',6,'FontName','Helvetica','location','n');
legend('boxoff');

set(gcf,'position',[0 0 350 150]);
set(gca,'position',[0.1 0.17 0.5 0.7]);