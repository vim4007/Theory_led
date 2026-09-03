
base_dir = fileparts(fileparts(fileparts(which('supplementary_fig1'))));
vir_tbl = readtable(fullfile(base_dir,'data','mouse','Scores','Virulence_screen_clean_table.csv'));
tbl     = readtable(fullfile(base_dir,'data','mouse','Scores','ProtectionScreen_CDI_mouse.csv'));
tbl = [tbl; vir_tbl];
tbl.cdiffstrain = string(tbl.cdiffstrain);
tbl(tbl.cdiffstrain=="vpi10463", :) = [];
tbl.experiment  = string(tbl.experiment);
tbl.exp_id      = tbl.experiment + "_" + tbl.cdiffstrain + "_" + string(tbl.mouse);

blue = [0 0.45 0.74];

%% 
primary = tbl(~contains(tbl.cdiffstrain,'.vpi') & ~contains(tbl.cdiffstrain,'.st1.75'), :);
secondary = tbl(endsWith(tbl.cdiffstrain,'.vpi') | strcmp(tbl.cdiffstrain,'vpi'), :);

%%
p = primary;
p.experiment  = categorical(p.experiment);
p.exp_id      = categorical(p.exp_id);
p.relweight(isnan(p.relweight)) = 0;
p.cdiffstrain = categorical(p.cdiffstrain);
p.cdiffstrain = reordercats(p.cdiffstrain, ['ui'; setdiff(categories(p.cdiffstrain),'ui')]);
p.relweight   = -1*p.relweight;  

vm = fitlme(p, 'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
ci = coefCI(vm);
V = table(vm.CoefficientNames(:), fixedEffects(vm), ci(:,1), ci(:,2), ...
    'VariableNames', {'name','Estimate','lo','hi'});
V(1,:) = [];  
V.Strain = string(extractAfter(V.name, "cdiffstrain_"));
V(strcmpi(V.Strain,'cd196'), :) = [];  
V = sortrows(V, 'Estimate', 'ascend');

figure('Color','w','Position',[100 100 1000 500]);
bar(V.Estimate, 'FaceColor', [0.2 0.45 0.8]); hold on
errorbar(1:height(V), V.Estimate, V.Estimate-V.lo, V.hi-V.Estimate, 'k.', 'LineWidth', 1.2);
set(gca,'XTick',1:height(V),'XTickLabel',upper(strrep(V.Strain,'.','-')),'XTickLabelRotation',45);
ylabel('Virulence Estimate'); ylim([-10 25]);
set(gca,'FontSize',25); box on; hold off

%%
p2 = primary; p2.cdiffstrain = string(p2.cdiffstrain);
strainsC = setdiff(unique(p2.cdiffstrain), ["", "cd196"]);   
strainsC = sort(strainsC);

figure('Color','w','Position',[60 60 1200 800]);
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
for s = 1:numel(strainsC)
    nexttile; hold on
    sdata = p2(p2.cdiffstrain==strainsC(s), :);
    [x,y,e,nTot,~] = strain_traj(sdata);
    if ~isempty(x)
        fill([x; flipud(x)],[y-e; flipud(y+e)], blue, 'FaceAlpha',0.3,'EdgeColor','none');
        plot(x, y, '-', 'Color', blue, 'LineWidth', 1.5);
    end
    title(sprintf('%s (%d mice)', upper(strrep(strainsC(s),'.','-')), nTot), 'FontSize', 9);
    ylim([60 120]); xlim([0 7]); box on; set(gca,'FontSize',15); hold off
end
title(tl, 'Primary infection weight trajectories');

%%
sec = secondary; sec.cdiffstrain = string(sec.cdiffstrain);
strainsD = setdiff(unique(sec.cdiffstrain), ["", "cd196.vpi"]);  
strainsD = sort(strainsD);

figure('Color','w','Position',[60 60 1200 800]);
tl = tiledlayout('flow','TileSpacing','compact','Padding','compact');
for s = 1:numel(strainsD)
    nexttile; hold on
    sdata = sec(sec.cdiffstrain==strainsD(s), :);
    [x,y,e,nTot,~] = strain_traj(sdata);
    if ~isempty(x)
        fill([x; flipud(x)],[y-e; flipud(y+e)], blue, 'FaceAlpha',0.3,'EdgeColor','none');
        plot(x, y, '-', 'Color', blue, 'LineWidth', 1.5);
    end
    label = upper(strrep(strrep(strainsD(s),'.vpi',''),'.','-'));
    title(sprintf('%s (%d mice)', label, nTot), 'FontSize', 9);
    ylim([60 120]); xlim([0 7]); box on; set(gca,'FontSize',15); hold off
end
title(tl, 'Co-infection (ST1 + VPI) weight trajectories');

%%
function [x,y,e,nTotal,nDied] = strain_traj(sdata)

    sdata.exp_id = string(sdata.exp_id);
    ids = unique(sdata.exp_id);

    isDied = false(numel(ids),1);
    for i = 1:numel(ids)
        md = sdata(sdata.exp_id==ids(i), :);
        isDied(i) = any(md.death==1) || any(md.relweight==0);
    end
    diedIds = ids(isDied);

    for i = 1:numel(diedIds)
        md = sortrows(sdata(sdata.exp_id==diedIds(i), :), 'day');
        last = find(md.relweight>0 & ~isnan(md.relweight), 1, 'last');
        if ~isempty(last)
            plot(md.day(1:last), md.relweight(1:last), 'k-', 'LineWidth', 1);
            plot(md.day(last), md.relweight(last), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
        end
    end

    days = sort(unique(sdata.day));
    y = nan(numel(days),1); e = nan(numel(days),1);
    for d = 1:numel(days)
        w = sdata.relweight(sdata.day==days(d));
        w = w(w>0 & ~isnan(w));
        if isempty(w), continue; end
        y(d) = mean(w);
        if numel(w) > 1, e(d) = std(w); end
    end
    e(isnan(e)) = 0;
    v = ~isnan(y); x = days(v); y = y(v); e = e(v);

    nTotal = numel(ids); nDied = numel(diedIds);
end