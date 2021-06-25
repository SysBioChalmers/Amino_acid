% load('GEM-yeast.mat');
% model_split = splitRevRxns(model);
% cd Models/;
% save('GEM-yeast-split.mat','model_split');
% cd ../;

load('GEM-yeast-split.mat');
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
enzymedata = updateYeastEnzyme(enzymedata,model_split);
enzymedata_org = enzymedata;

step = 0.001;
mu = 1;
tol = mu*0.00001;
substrateExID = 'r_1714'; % glucose

% max sa for estimating cost
model_max = convertModel(model_split,enzymedata);
model_max = setYeastMedia(model_max);
% model_max = blockYeastRxns(model_max);
cost_maxSA = calculateYeastCost(model_max,substrateExID,mu,step,tol);

cd YeastCost/;
save('cost_maxSA.mat','cost_maxSA');
cd ../;


% condition-specific sa for estimating cost
load('YeastSA.mat');

cond = YeastSA.condition(~contains(YeastSA.condition,'etoh','IgnoreCase',true));

prot_cost = zeros(20,length(cond));

for k = 1:length(cond)
    display([num2str(k),'/',num2str(length(cond))]);
    enzymedata = enzymedata_org;
    enzymedata = updateCondSpecYeastEnzyme(enzymedata,cond(k));
    
    model = convertModel(model_split,enzymedata);
    
    model = setYeastMedia(model);
    model = blockYeastRxns(model);
    
    cost_CondSA = calculateYeastCost(model,substrateExID,mu,step,tol);
    
    cd YeastCost/;
    save(['cost_CondSA_' cond{k} '.mat'],'cost_CondSA');
    cd ../;
    
    
    
%     flux_ref = cost_yeast.flux_ref;
%     flux_AAs = cost_yeast.flux_AAs;
%     
%     tot_enzymatic_rxns = zeros(length(AA),1);
%     
%     figure();
%     for i = 1:length(AA)
%         change = flux_AAs(:,i)-flux_ref;
%         change(abs(change)<1e-7) = 0;
%         active_rxn = model.rxns(change ~= 0);
%         active_flux = change(change ~= 0);
%         
%         enzymatic_rxn = intersect(active_rxn,enzymedata.rxn);
%         [~,b] = ismember(enzymatic_rxn,enzymedata.rxn);
%         enzymatic_kcat_conf = enzymedata.kcat_conf(b);
%         %     enzymatic_kcat = enzymedata.kcat(b);
%         %     printRxnFormula(model,'rxnAbbrList',enzymatic_rxn,'metNameFlag',true);
%         x = [0 1 2 3 4 5 6];
%         y = zeros(1,7);
%         for j = 1:length(x)
%             y(j) = sum(enzymatic_kcat_conf == x(j));
%         end
%         pert = round(sum(y(x>5))/sum(y)*100,0);
%         pert_text = [num2str(pert) '%'];
%         tot_enzymatic_rxns(i) = sum(y);
%         
%         subplot(4,5,i);
%         b = bar(x,y,0.7,'FaceColor','flat','LineWidth',0.5);
%         b.CData(1:6,:) = repmat([64,64,64]/255,6,1);
%         b.CData(7,:) = repmat([202,0,32]/255,1,1);
%         b.EdgeColor = 'w';
%         set(gca,'XColor','k');
%         set(gca,'YColor','k');
%         set(gca,'FontSize',6,'FontName','Helvetica');
%         
%         if max(y) < 25
%             ylim([0 30]);
%             text(2,max(y),pert_text,'FontSize',7,'FontName','Helvetica','Color',[202,0,32]/255);
%             text(-0.5,27,AA{i},'FontSize',7,'FontName','Helvetica','Color','k','FontWeight','bold');
%         else
%             ylim([0 60]);
%             text(2,max(y),pert_text,'FontSize',7,'FontName','Helvetica','Color',[202,0,32]/255);
%             text(-0.5,54,AA{i},'FontSize',7,'FontName','Helvetica','Color','k','FontWeight','bold');
%         end
%         ylabel('count','FontSize',7,'FontName','Helvetica','Color','k');
%         xlabel('confidence score','FontSize',7,'FontName','Helvetica','Color','k');
%     end
%     set(gcf,'position',[200 200 450 350]);
    
    
    
    
    
end


