
base_dir = fileparts(fileparts(fileparts(which('supp_fig3'))));
tbl = readtable(fullfile(base_dir, 'data', '16s_sequencing', 'tblAbund.xls'));
tbl.Initial_infection   = string(tbl.Initial_infection);
tbl.Initial_antibiotics = string(tbl.Initial_antibiotics);

genusCols = tbl.Properties.VariableNames(6:end);

meanAbund = mean(tbl{:, genusCols}, 1);
[~, sortIdx] = sort(meanAbund, 'descend');
namedIdx = sortIdx(~strcmp(genusCols(sortIdx), 'Other'));
top9     = genusCols(namedIdx(1:9));
otherIdx = namedIdx(10:end);
tbl.Other_combined = tbl.Other + sum(tbl{:, genusCols(otherIdx)}, 2);
plotCols   = [top9, {'Other_combined'}];
plotLabels = [top9, {'Other'}];
nPlot      = numel(plotCols);

groupKeys   = {'uninfected', 'ST1_75', 'ST1_12'};
groupLabels = {'Uninfected', 'ST1.75', 'ST1.12'};
groupColors = [0.5 0.5 0.5; 0 0.45 0.74; 0.9 0.5 0.3];
days = [1 2 7];

hexColors = {'#5875DE','#1CF8EC','#2A9D8F','#0D7E2B','#CA0BE8', ...
             '#FBA22E','#BEA89A','#7D6E65','#1A7A73','#D4B896'};
genusPalette = cell2mat(cellfun(@(h) sscanf(h(2:end),'%2x%2x%2x')'/255, ...
                        hexColors(:), 'UniformOutput', false));

%% 
figure('Color','w','Position',[80 80 1150 450]);
tlA = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for d = 1:numel(days)
    nexttile; 
    sub = tbl(tbl.Day==days(d) & tbl.Initial_antibiotics=="mnvc", :);
    barData = zeros(3, nPlot);
    for g = 1:3
        grp = sub(sub.Initial_infection==groupKeys{g}, :);
        if ~isempty(grp), barData(g,:) = mean(grp{:, plotCols},1)/100; end
    end
    b = bar(barData, 'stacked', 'EdgeColor','none');
    for i = 1:nPlot, b(i).FaceColor = genusPalette(i,:); end
    set(gca,'XTick',1:3,'XTickLabel',groupLabels,'XTickLabelRotation',30,'FontSize',20);
    ylim([0 1]); title(sprintf('Day %d', days(d)));
    if d==1, ylabel('Relative abundance'); end
    box on;
end
lgd = legend(plotLabels, 'Location','eastoutside'); set(lgd,'Interpreter','none','FontSize',12);
lgd.Layout.Tile = 'east';
title(tlA, 'Relative abundance of gut taxa (mnvc)');

%% 
figure('Color','w','Position',[80 80 1150 400]);
tlB = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for d = 1:numel(days)
    nexttile; hold on
    sub = tbl(tbl.Day==days(d) & tbl.Initial_antibiotics=="mnvc" & ...
             (tbl.Initial_infection=="ST1_75" | tbl.Initial_infection=="ST1_12"), :);
    abd = sub{:, genusCols}/100;
    nS = height(sub);

    BC = zeros(nS);
    for i = 1:nS
        for j = i+1:nS
            v = sum(abs(abd(i,:)-abd(j,:))) / (sum(abd(i,:))+sum(abd(j,:)));
            BC(i,j)=v; BC(j,i)=v;
        end
    end
    [coords, eig] = cmdscale(BC);
    varExp = eig / sum(eig(eig>0)) * 100;

    i75 = sub.Initial_infection=="ST1_75";
    i12 = sub.Initial_infection=="ST1_12";
    scatter(coords(i75,1), coords(i75,2), 120, groupColors(2,:), 'filled', ...
        'MarkerFaceAlpha',0.8, 'DisplayName','ST1.75');
    scatter(coords(i12,1), coords(i12,2), 120, groupColors(3,:), 'filled', ...
        'MarkerFaceAlpha',0.8, 'DisplayName','ST1.12');
    xlabel(sprintf('PCo1 (%.1f%%)', varExp(1)));
    ylabel(sprintf('PCo2 (%.1f%%)', varExp(2)));
    title(sprintf('PCoA Bray-Curtis - Day %d', days(d)));
    set(gca,'FontSize',20); box on; hold off
end
lgd = legend({'ST1.75','ST1.12'}, 'Location','eastoutside');
lgd.Layout.Tile = 'east';
title(tlB, 'Beta diversity (Bray-Curtis PCoA)');