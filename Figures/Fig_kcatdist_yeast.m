expID = 'Yu_Clim';
% expID = 'DiBartolomeo_Gluc';

load(['cost_yeast_' expID '.mat']);
load(['enzymedataYeast_' expID '.mat']);
load(['modelYeast_' expID '.mat']);
aaList = {'ala';'arg';'asn';'asp';'cys';
'gln';'glu';'gly';'his';'ile';
'leu';'lys';'met';'phe';'pro';
'ser';'thr';'trp';'tyr';'val'};
AA = cost_yeast.AA;
flux_ref = cost_yeast.flux_ref;
flux_AAs = cost_yeast.flux_AAs;

tot_enzymatic_rxns = zeros(length(AA),1);

figure();
for i = 1:length(AA)
    change = flux_AAs(:,i)-flux_ref;
    change(abs(change)<1e-7) = 0;
    active_rxn = model.rxns(change ~= 0);
    active_flux = change(change ~= 0);
    
    enzymatic_rxn = intersect(active_rxn,enzymedata.rxn);
    [~,b] = ismember(enzymatic_rxn,enzymedata.rxn);
    enzymatic_kcat_conf = enzymedata.kcat_conf(b);
%     enzymatic_kcat = enzymedata.kcat(b);
%     printRxnFormula(model,'rxnAbbrList',enzymatic_rxn,'metNameFlag',true);
    x = [0 1 2 3 4 5];
    y = zeros(1,6);
    for j = 1:length(x)
        y(j) = sum(enzymatic_kcat_conf == x(j));
    end
    pert = round(sum(y(x>4))/sum(y)*100,0);
    pert_text = [num2str(pert) '%'];
    tot_enzymatic_rxns(i) = sum(y);
    
    subplot(4,5,i);
    b = bar(x,y,0.7,'FaceColor','flat','LineWidth',0.5);
    b.CData(1:5,:) = repmat([64,64,64]/255,5,1);
    b.CData(6,:) = repmat([202,0,32]/255,1,1);
    b.EdgeColor = 'w';
    set(gca,'XColor','k');
    set(gca,'YColor','k');
    set(gca,'FontSize',6,'FontName','Helvetica');
    
    if max(y) < 25
        ylim([0 30]);
        text(2,max(y),pert_text,'FontSize',7,'FontName','Helvetica','Color',[202,0,32]/255);
        text(-0.5,27,aaList{i},'FontSize',7,'FontName','Helvetica','Color','k','FontWeight','bold');
    else
        ylim([0 60]);
        text(2,max(y),pert_text,'FontSize',7,'FontName','Helvetica','Color',[202,0,32]/255);
        text(-0.5,54,aaList{i},'FontSize',7,'FontName','Helvetica','Color','k','FontWeight','bold');
    end
    ylabel('count','FontSize',7,'FontName','Helvetica','Color','k');
    xlabel('confidence score','FontSize',7,'FontName','Helvetica','Color','k');
end
set(gcf,'position',[200 200 450 350]);

