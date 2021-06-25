load('GEM-yeast-split.mat');
rxnlist = model_split.rxns;
grlist = model_split.grRules;
rowlist = cell(0,1);
salist = zeros(length(rxnlist),0);
file = dir('*.mat');
for i = 1:length(file)
    filename = file(i).name;
    if isempty(regexp(filename,'R\d*.mat','match'))
        load(filename);
        if YeastSA.deviation < 0.15
            condid = strrep(filename,'YeastSA_','');
            condid = strrep(condid,'.mat','');
            rowlist = [rowlist {condid}];
            satmp = zeros(length(rxnlist),1);
            [~,b] = ismember(YeastSA.rxn,rxnlist);
            satmp(b) = YeastSA.sa;
            salist = [salist satmp];
        end
    end
end


idx1 = ~any(salist == 0,2); % remain rxns that have sa at any conditions
% idx1 = any(salist,2); % remain rxns that have sa at least one condition
rxnlist = rxnlist(idx1);
grlist = grlist(idx1);
salist = salist(idx1,:);

% calculate kapp
idx2 = ~contains(grlist,'and');
rxnlist_kapp = rxnlist(idx2);
grlist_kapp = grlist(idx2);
salist_kapp = salist(idx2,:);

[num,txt,~] = xlsread('UniProt_Yeast.xlsx');
geneUniProt = txt(2:end,1);
MWUniProt = num;
clear num txt;

kapplist = zeros(length(grlist_kapp),length(rowlist));
for i = 1:length(grlist_kapp)
    
    gr_tmp = grlist_kapp{i};
    gr_tmp = strrep(gr_tmp,'(','');
    gr_tmp = strrep(gr_tmp,')','');
    gr_tmp = strrep(gr_tmp,' or ',' ');
    gr_tmp = strtrim(gr_tmp);
    gr_tmp = strsplit(gr_tmp);
    gr_tmp = unique(gr_tmp)';
    
    [~,idxtmp] = ismember(gr_tmp,geneUniProt);
    
    kapplist(i,:) = salist_kapp(i,:) * max(MWUniProt(idxtmp)) / 1000 / 60;
    
end

[m,n] = size(kapplist);

[RHO,~] = corr(log10(kapplist),'type','Pearson');
figure();
xlbl = num2cell(1:1:length(rowlist));
xlbl_tmp = cellfun(@(x) [num2str(x),'.'],xlbl,'UniformOutput',false);
ylbl = strcat(xlbl_tmp,rowlist);
minclr = [255,255,178]/255;
maxclr = [242,94,13]/255;
tmp1 = linspace(minclr(1),maxclr(1),129)';
tmp2 = linspace(minclr(2),maxclr(2),129)';
tmp3 = linspace(minclr(3),maxclr(3),129)';
clrmap = [tmp1 tmp2 tmp3];
h = heatmap(xlbl,ylbl,RHO,'Colormap',clrmap,...
    'ColorMethod','count','CellLabelColor','k');
title(['Correlations between kapps (N = ' num2str(m) ') under ' num2str(n) ' conditions']);
set(h,'FontSize',6,'FontName','Helvetica');
h.FontColor = 'k';
set(gcf,'position',[320 320 260 260]);
set(gca,'position',[0.23 0.2 0.6 0.6]);

% correlate kmax to kcat
enzymedata = collectkcats(model_split.rxns(~ismember(model_split.rules,{''})),model_split,'saccharomyces cerevisiae');
idx3 = enzymedata.kcat_conf > 2;
rxnlist_kcat = enzymedata.rxn(idx3);
kcatlist = enzymedata.kcat(idx3)/3600;
kmaxlist = max(kapplist,[],2);
rxns = intersect(rxnlist_kcat,rxnlist_kapp);

[~,p] = ismember(rxns,rxnlist_kcat);
kcat = kcatlist(p);
[~,q] = ismember(rxns,rxnlist_kapp);
kmax = kmaxlist(q);

figure();
[RHO,PVAL] = corr(log10(kcat),log10(kmax),'Type','Pearson');
scatter(log10(kcat),log10(kmax),10,'o','filled','LineWidth',1,'MarkerEdgeColor',[1,1,1],'MarkerFaceColor',[0,0,0],'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0);
xlim([-2 6]);
ylim([-2 6]);
text(3.5,-1.5,['N' ' = ' num2str(length(kcat))],'FontSize',7,'FontName','Helvetica');
text(-1.5,5.7,['R^2' '= ' num2str(round(RHO^2,3))],'FontSize',7,'FontName','Helvetica');
text(-1.5,4.7,['p' ' = ' num2str(round(PVAL,2))],'FontSize',7,'FontName','Helvetica');
title('S. cerevisiae','FontSize',7,'FontName','Helvetica');
set(gca,'FontSize',6,'FontName','Helvetica');
xlabel('log10 kcat (in vitro) [/s]','FontSize',7,'FontName','Helvetica');
ylabel('log10 kmax (in vivo) [/s]','FontSize',7,'FontName','Helvetica');
set(gcf,'position',[200 200 120 120]);
set(gca,'position',[0.2 0.2 0.75 0.75]);
