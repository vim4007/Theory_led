
try
    thisFile = mfilename('fullpath');
    if isempty(thisFile), error('run-as-section'); end
catch
    thisFile = matlab.desktop.editor.getActiveFilename;   % running a section
end
base_dir = fileparts(fileparts(fileparts(thisFile)));

UI_mice        = [1 6 7];
AVIR_mice      = [11 12 13];
GRP_COLORS     = [0.4 0.6 0.9; 0.9 0.5 0.3];   % UI, ST1-75 colonized

%% 
A_types  = {'B_cells','CD19negTCRbneg','NK_cells','T_cells','CD4pos','CD8pos'};
A_pretty = {'B cells','CD19^{-}TCR\beta^{-}','NK cells','T cells','CD4^{+}','CD8^{+}'};
A_frac = flow_fractions(base_dir, 'ks10_adaptive_csv_files', ...
                        'all_adaptive/all_adaptive', A_types, 14);
plot_immune_bars(A_frac, A_types, A_pretty, UI_mice, AVIR_mice, GRP_COLORS, ...
    'Adaptive immune cell type', 'Adaptive immune fractions');

%% 
B_types  = {'DCs','cd11b_pos','lymphoid_DCs','macrophages','monocytes','myeloid_DCs','neutrophils'};
B_pretty = {'DCs','CD11b^{+}','Lymphoid DCs','Macrophages','Monocytes','Myeloid DCs','Neutrophils'};
B_frac = flow_fractions(base_dir, 'ks10_innate_csv_files', ...
                        'all_cells/all_cells', B_types, 13);
plot_immune_bars(B_frac, B_types, B_pretty, UI_mice, AVIR_mice, GRP_COLORS, ...
    'Innate immune cell type', 'Innate immune fractions');

%% 
C = readtable(fullfile(base_dir,'data','mouse','rag1ko','KS11_rechallenge_relweight.xlsx'));
C.group = string(C.group); C.mouse_id = string(C.mouse_id);

cKeys   = {'uninfected_b6','st175_b6','st175_rag1ko'};
cTitles = {'Uninfected B6','ST1-75 B6','ST1-75 RAG1 KO'};
cColors = [0.50 0.50 0.50; 0.00 0.45 0.85; 0.80 0.10 0.55];

figure('Color','w','Position',[100 100 640 460]); hold on
for g = 1:3
    grp = C(C.group == cKeys{g}, :);
    ids = unique(grp.mouse_id);
    surv = ids(arrayfun(@(m) ~any(grp.relweight(grp.mouse_id==m & grp.day>0)==0 | ...
                        isnan(grp.relweight(grp.mouse_id==m & grp.day>0))), ids));
    days = sort(unique(grp.day));
    mu = nan(numel(days),1); se = nan(numel(days),1);
    for d = 1:numel(days)
        vals = grp.relweight(ismember(grp.mouse_id,surv) & grp.day==days(d));
        vals = vals(~isnan(vals));
        if ~isempty(vals)
            mu(d) = mean(vals);
            if numel(vals) > 1, se(d) = std(vals)/sqrt(numel(vals)); end
        end
    end
    se(isnan(se)) = 0; v = ~isnan(mu);
    x = days(v); y = mu(v); e = se(v);
    fill([x; flipud(x)], [y-e; flipud(y+e)], cColors(g,:), ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    plot(x, y, '-o', 'Color',cColors(g,:), 'MarkerFaceColor',cColors(g,:), ...
        'LineWidth',3.5, 'MarkerSize',5, 'DisplayName',cTitles{g});
    fprintf('Group: %-16s | Surviving: %d\n', cTitles{g}, numel(surv));
end
yline(100,'k:','HandleVisibility','off');
xlabel('Day after secondary challenge'); ylabel('Relative weight (%)');
legend('Location','southwest','Box','off','FontSize',20);
xlim([0 7]); ylim([70 115]); set(gca,'FontSize',20); box on; hold off
title('Protection in RAG1 knockout mice','FontSize',15);

%% 
M = readtable(fullfile(base_dir,'data','16s_sequencing','tblAbund.xls'));
M.Initial_infection   = string(M.Initial_infection);
M.Initial_antibiotics = string(M.Initial_antibiotics);

genusCols = M.Properties.VariableNames(6:end);
[~, sIdx]  = sort(mean(M{:,genusCols},1),'descend');
top12      = genusCols(sIdx(1:12));
resid_genus = genusCols(~strcmpi(genusCols,'Clostridioides'));

shannon = @(p) -sum(p(p>0).*log(p(p>0)));
M.ResidShannon = nan(height(M),1);
for i = 1:height(M)
    v = M{i,resid_genus}; s = sum(v);
    if s > 0, M.ResidShannon(i) = shannon(v/s); end
end

days = [0 1 2 7];
blue = [0 0.45 0.85]; orange = [0.90 0.50 0.30];
mean75=nan(1,4); sem75=nan(1,4); mean12=nan(1,4); sem12=nan(1,4); p_shan=nan(1,4);

for d = 1:4
    sub = M(M.Day==days(d) & M.Initial_antibiotics=='mnvc', :);
    s75 = sub.ResidShannon(sub.Initial_infection=='ST1_75'); s75 = s75(~isnan(s75));
    s12 = sub.ResidShannon(sub.Initial_infection=='ST1_12'); s12 = s12(~isnan(s12));
    if ~isempty(s75), mean75(d)=mean(s75); sem75(d)=std(s75)/sqrt(numel(s75)); end
    if ~isempty(s12), mean12(d)=mean(s12); sem12(d)=std(s12)/sqrt(numel(s12)); end
    if ~isempty(s75)&&~isempty(s12), p_shan(d)=ranksum(s75,s12); end
end
q_shan = bh_fdr(p_shan);

p_bc = nan(1,4); rng(42); nPerm = 5000;
for d = 1:4
    sub = M(M.Day==days(d) & M.Initial_antibiotics=='mnvc' & ...
           (M.Initial_infection=='ST1_75'|M.Initial_infection=='ST1_12'), :);
    if height(sub) < 4, continue; end
    Ab = sub{:,resid_genus}; Ab = Ab./sum(Ab,2);
    grp = double(sub.Initial_infection=='ST1_75'); nS = size(Ab,1);
    BC = zeros(nS);
    for i=1:nS, for j=i+1:nS
        d_ij = sum(abs(Ab(i,:)-Ab(j,:)))/(sum(Ab(i,:))+sum(Ab(j,:)));
        BC(i,j)=d_ij; BC(j,i)=d_ij;
    end, end
    stat = @(g) mean(BC((g~=g'))) - mean(BC((g==g')&~eye(nS)));
    obs = stat(grp);
    sp = arrayfun(@(~) stat(grp(randperm(nS))), 1:nPerm);
    p_bc(d) = (sum(sp>=obs)+1)/(nPerm+1);
end
q_bc = bh_fdr(p_bc);

sub1 = M(M.Day==1 & M.Initial_antibiotics=='mnvc' & ...
        (M.Initial_infection=='ST1_75'|M.Initial_infection=='ST1_12'), :);
Aall = sub1{:,resid_genus}; Aall = Aall./sum(Aall,2);
rtop = top12(ismember(top12,resid_genus));
[~,cpos] = ismember(rtop,resid_genus);
g75 = sub1.Initial_infection=='ST1_75'; g12 = sub1.Initial_infection=='ST1_12';
p_gen = nan(1,numel(rtop));
for k = 1:numel(rtop)
    p_gen(k) = ranksum(Aall(g75,cpos(k)), Aall(g12,cpos(k)));
end
q_gen = bh_fdr(p_gen);
min_q = min([q_shan q_bc q_gen],[],'omitnan');

figure('Color','w','Position',[100 100 720 500]); hold on
errorbar(days, mean75, sem75, '-o','Color',blue,'MarkerFaceColor',blue, ...
    'LineWidth',2.5,'MarkerSize',7,'CapSize',6,'DisplayName','ST1-75');
errorbar(days, mean12, sem12, '-o','Color',orange,'MarkerFaceColor',orange, ...
    'LineWidth',2.5,'MarkerSize',7,'CapSize',6,'DisplayName','ST1-12');
yl = ylim;
for d = 1:4
    if ~isnan(q_shan(d))
        ytop = max(mean75(d)+sem75(d), mean12(d)+sem12(d));
        text(days(d), ytop+0.06*range(yl), sprintf('q=%.2f',q_shan(d)), ...
            'HorizontalAlignment','center','FontSize',11,'Color',[0.3 0.3 0.3]);
    end
end
txt = sprintf('Day 1 Bray-Curtis p=%.3f\nTop-genus min q=%.2f', p_bc(2), min(q_gen(~isnan(q_gen))));
text(xlim*[0.97;0.03], yl(1)+0.10*range(yl), txt, 'FontSize',15,'VerticalAlignment','bottom');
xlabel('Day'); ylabel('Residual Shannon index');
legend('Location','northwest','Box','off','FontSize',20);
xlim([-0.3 7.3]); set(gca,'XTick',days,'FontSize',20); box on; hold off
title('Residual microbiota composition','FontSize',15);

fprintf('\n--- Shannon (per day, BH) ---\n');
for d=1:4, fprintf('  Day %d: p=%.3f q=%.3f\n', days(d), p_shan(d), q_shan(d)); end
fprintf('--- Bray-Curtis (per day, BH) ---\n');
for d=1:4, fprintf('  Day %d: p=%.3f q=%.3f\n', days(d), p_bc(d), q_bc(d)); end
fprintf('--- Top-12 genera at day 1 (BH) ---\n');
for k=1:numel(rtop), fprintf('  %-20s p=%.3f q=%.3f\n', rtop{k}, p_gen(k), q_gen(k)); end
fprintf('\nOverall minimum q = %.2f\n', min_q);

%% 
G = readtable(fullfile(base_dir,'data','Genomics','binary_gene_pre_abs.csv'));
strainNames = string(G{:,1});
geneData    = double(G{:,2:end-2});
protection  = G{:,end};

accessory = geneData(:, ~all(geneData==1,1));       % drop core genes
[~, score, ~, ~, explained] = pca(accessory - mean(accessory,1));

protCmap = [zeros(256,1), linspace(0.8,0,256)', linspace(0,1,256)'];
minP = min(protection); maxP = max(protection);

figure('Color','w','Position',[100 100 600 500]); hold on
for i = 1:numel(strainNames)
    ci = round((protection(i)-minP)/(maxP-minP)*255)+1;
    scatter(score(i,1), score(i,2), 250, protCmap(ci,:), 'filled', ...
        'MarkerEdgeColor','k','LineWidth',0.5);
    text(score(i,1)+0.3, score(i,2)+0.3, strrep(strainNames{i},'ST1-','ST1.'), ...
        'FontSize',25,'Interpreter','none');
end
colormap(protCmap); cb = colorbar; cb.Label.String = 'Protection score';
caxis([minP maxP]);
xlabel(sprintf('PC1 (%.1f%%)',explained(1)));
ylabel(sprintf('PC2 (%.1f%%)',explained(2)));
title('Accessory-gene PCA','FontSize',15);
set(gca,'FontName','Arial','FontSize',20); box on; hold off

%% ============================================================
%  Local functions
%  ============================================================
function frac = flow_fractions(base_dir, csvRoot, allCellsStub, cellTypes, nMetrics)
    data_dir = fullfile(base_dir, 'data', 'flow_cytometry');
    mice = [1 6 7 11 12 13];

    allMice = [];
    for i = 1:numel(mice)
        t = readtable(fullfile(data_dir, sprintf('%s/%s_Specimen_001_%d.csv', ...
            csvRoot, allCellsStub, mice(i))));
        t.mouse(:) = mice(i);
        allMice = [allMice; t];
    end
    eventVars = setdiff(allMice.Properties.VariableNames, {'mouse'}, 'stable');

    gated = allMice;
    for j = 1:numel(cellTypes)
        d = dir(fullfile(data_dir, sprintf('%s/%s/*.csv', csvRoot, cellTypes{j})));
        oftype = [];
        for i = 1:numel(d)
            t = readtable(fullfile(data_dir, sprintf('%s/%s/%s', csvRoot, cellTypes{j}, d(i).name)));
            oftype = [oftype; t];
        end
        gated = addvars(gated, ismember(allMice(:,eventVars), oftype(:,eventVars)));
        gated.Properties.VariableNames{end} = cellTypes{j};
    end
    gCols = (width(gated)-numel(cellTypes)+1) : width(gated);
    passed = gated(sum(gated{:,gCols},2) > 0, :);
    lo = min(passed{:,1:nMetrics}); hi = max(passed{:,1:nMetrics});
    inRange = all(gated{:,1:nMetrics} >= lo & gated{:,1:nMetrics} <= hi, 2);
    cells = gated(inRange, :);

    ids = unique(cells.mouse);
    frac = table(ids, 'VariableNames', {'mouse'});
    for i = 1:numel(ids)
        dm = cells(cells.mouse==ids(i), :);
        frac{i, 2:numel(cellTypes)+1} = sum(dm{:,cellTypes},1) / height(dm);
    end
    frac.Properties.VariableNames(2:end) = cellTypes;
end

function plot_immune_bars(frac, cellTypes, pretty, UI, AVIR, colors, xlab, ttl)

    ids = frac.mouse; vals = frac{:,2:end};
    UI_d = vals(ismember(ids,UI), :);  AV_d = vals(ismember(ids,AVIR), :);
    means = [mean(UI_d,1); mean(AV_d,1)];
    sds   = [std(UI_d,0,1); std(AV_d,0,1)];
    n = numel(cellTypes);

    p = arrayfun(@(j) ranksum(UI_d(:,j), AV_d(:,j)), 1:n);
    fprintf('\n%s Wilcoxon p-values (n=3/group):\n', ttl);
    for j=1:n, fprintf('  %s: p=%.4f\n', cellTypes{j}, p(j)); end

    gw = 0.7; bw = gw/2;
    figure('Color','w'); hold on
    for g = 1:2
        xp = (1:n) - gw/2 + (g-0.5)*bw;
        bar(xp, means(g,:), bw, 'FaceColor', colors(g,:), 'EdgeColor','none','FaceAlpha',0.85);
    end
    for g = 1:2
        xp = (1:n) - gw/2 + (g-0.5)*bw;
        gd = UI_d; if g==2, gd = AV_d; end
        errorbar(xp, means(g,:), sds(g,:), 'k.', 'LineWidth',1.5, 'CapSize',5);
        for j = 1:n
            scatter(repmat(xp(j),size(gd,1),1), gd(:,j), 40, 'k', 'filled', 'MarkerFaceAlpha',0.6);
        end
    end
    for j = 1:n
        x1 = j-gw/2+0.5*bw; x2 = j-gw/2+1.5*bw;
        ymax = max(means(:,j)) + max(sds(:,j)) + 0.005;
        star = ''; if p(j)<0.001, star='***'; elseif p(j)<0.01, star='**'; elseif p(j)<0.05, star='*'; end
        if ~isempty(star)
            plot([x1 x2],[ymax ymax],'k-','LineWidth',1);
            text((x1+x2)/2, ymax+0.002, star, 'FontSize',16,'FontWeight','bold','HorizontalAlignment','center');
        end
    end
    set(gca,'XTick',1:n,'XTickLabel',pretty,'XTickLabelRotation',45,'TickLabelInterpreter','tex');
    legend({'Uninfected (UI)','ST1-75 colonized'},'Location','northeast');
    set(gca,'FontName','Arial','FontSize',20);
    xlabel(xlab); ylabel('Fraction of gated cells');
    title(ttl,'FontSize',15); box on; grid on; hold off
end

function q = bh_fdr(p)
    q = nan(size(p));
    valid = ~isnan(p);
    pv = p(valid); [ps, si] = sort(pv(:)'); m = numel(ps);
    if m == 0, return; end
    qs = ps .* m ./ (1:m);
    qs = min(1, fliplr(cummin(fliplr(qs))));
    qt = nan(1,m); qt(si) = qs;
    q(valid) = qt;
end