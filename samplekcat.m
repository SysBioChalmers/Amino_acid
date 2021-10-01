
load('GEM-yeast-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);
model = convertModel(model_split,enzymedata);
clear model_split;

step = 0.0001;
mu = 0.4;

n_sample = 1000;
idx_changedkcat = enzymedata.kcat_conf < 5;
n_changedkcat = sum(idx_changedkcat);
exponent = 2;
digit = 0.001;
random_data = 10.^[randi([-exponent/digit,exponent/digit],n_changedkcat,n_sample)*digit];
% random_data = ones(n_changedkcat,n_sample);

load('YeastCost/cost_kmax.mat');
cost_random.AA = cost_kmax.AA;
cost_random.cost_prot = zeros(20,n_sample);

for k = 1:n_sample
    ezmdata_rxn = enzymedata.rxn;
    ezmdata_minPC = enzymedata.minPC;
    ezmdata_minPC(idx_changedkcat) = ezmdata_minPC(idx_changedkcat).*random_data(:,k);
    
    disp([num2str(k),'/',num2str(n_sample)]);
    
    for i = 1:length(cost_random.AA)
        
        fluxAA = cost_kmax.flux_AAs(:,i);
        
        flux_tmp = fluxAA - cost_kmax.flux_ref;
        active_rxns = model.rxns(flux_tmp ~= 0);
        active_fluxes = flux_tmp(flux_tmp ~= 0);
        protein_cost = zeros(length(active_rxns),1);
        
        for j = 1:length(active_rxns)
            if ismember(active_rxns(j),ezmdata_rxn)
                protein_cost(j) = ezmdata_minPC(ismember(ezmdata_rxn,active_rxns(j)));
            end
        end
        cost_random.cost_prot(i,k) = sum(active_fluxes.*protein_cost)/step/mu/1000;
    end
    
end

cost_random.r = zeros(1,n_sample);
load('YeastAA.mat');
for i = 1:n_sample
    [RHO,~] = corr(cost_random.cost_prot(:,i),log(YeastAA.AAcontent_avg*100),'Type','Pearson');
    cost_random.r(i) = RHO;
end

[f,xi] = ksdensity(cost_random.r);
hold on;
plot(xi,f);
line([-0.5,-0.5],[0,max(f)],'Color','red');
text(-0.45,max(f),[num2str(100*sum(cost_random.r>-0.5)/n_sample) '% on the right'],'VerticalAlignment','middle','Color','red','FontSize',6,'FontName','Helvetica');
xlim([-1 ceil(10*max(cost_random.r))/10]);
title(['Distribution of ' num2str(n_sample) ' r values'],'FontSize',6,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('Pearson r','FontSize',7,'FontName','Helvetica');
ylabel('Density','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[600 400 150 100]);
set(gca,'position',[0.27 0.27 0.62 0.62]);

T = table(transpose(cost_random.r));
writetable(T,'1000_r_values.txt');


