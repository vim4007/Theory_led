base_dir = fileparts(fileparts(fileparts(which('supp_fig2_tsne'))));

UI_mice   = [1 6 7];
AVIR_mice = [11 12 13];

adaptive.root    = 'ks10_adaptive_csv_files';
adaptive.stub    = 'all_adaptive/all_adaptive';
adaptive.types   = {'B_cells','CD19negTCRbneg','NK_cells','T_cells','CD4pos','CD8pos'};
adaptive.pretty  = {'B cells','CD19^{-}TCR\beta^{-}','NK cells','T cells','CD4^{+}','CD8^{+}'};
adaptive.nMetric = 14;

innate.root    = 'ks10_innate_csv_files';
innate.stub    = 'all_cells/all_cells';
innate.types   = {'DCs','cd11b_pos','lymphoid_DCs','macrophages','monocytes','myeloid_DCs','neutrophils'};
innate.pretty  = {'DCs','CD11b^{+}','Lymphoid DCs','Macrophages','Monocytes','Myeloid DCs','Neutrophils'};
innate.nMetric = 13;

Ad = tsne_panel(base_dir, adaptive, UI_mice, AVIR_mice);
In = tsne_panel(base_dir, innate,   UI_mice, AVIR_mice);

%%
figure('Color','w','Position',[80 80 1300 1000]);
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile; plot_by_celltype(Ad); title('Adaptive: cell types');
nexttile; plot_by_group(Ad);    title('Adaptive: infection group');
nexttile; plot_by_celltype(In); title('Innate: cell types');
nexttile; plot_by_group(In);    title('Innate: infection group');

title(tl, 'tSNE of flow-cytometry profiles','FontWeight','bold');

%%
function P = tsne_panel(base_dir, cfg, UI_mice, AVIR_mice)
    data_dir = fullfile(base_dir, 'data', 'flow_cytometry');
    mice = [1 6 7 11 12 13];

    allMice = [];
    for i = 1:numel(mice)
        t = readtable(fullfile(data_dir, sprintf('%s/%s_Specimen_001_%d.csv', ...
            cfg.root, cfg.stub, mice(i))));
        t.mouse(:) = mice(i);
        allMice = [allMice; t];
    end
    eventVars = setdiff(allMice.Properties.VariableNames, {'mouse'}, 'stable');

    gated = allMice;
    for j = 1:numel(cfg.types)
        d = dir(fullfile(data_dir, sprintf('%s/%s/*.csv', cfg.root, cfg.types{j})));
        oftype = [];
        for i = 1:numel(d)
            t = readtable(fullfile(data_dir, sprintf('%s/%s/%s', cfg.root, cfg.types{j}, d(i).name)));
            oftype = [oftype; t];
        end
        gated = addvars(gated, ismember(allMice(:,eventVars), oftype(:,eventVars)));
        gated.Properties.VariableNames{end} = cfg.types{j};
    end

    gCols = (width(gated)-numel(cfg.types)+1) : width(gated);
    passed = gated(sum(gated{:,gCols},2) > 0, :);
    lo = min(passed{:,1:cfg.nMetric}); hi = max(passed{:,1:cfg.nMetric});
    inRange = all(gated{:,1:cfg.nMetric} >= lo & gated{:,1:cfg.nMetric} <= hi, 2);
    cells = gated(inRange, :);

    lab = repmat("Other", height(cells), 1);
    for j = 1:numel(cfg.types)
        lab(cells{:,cfg.types{j}}==1) = string(cfg.types{j});
    end
    P.cellLabels = categorical(lab);

    grp = repmat("Other", height(cells), 1);
    grp(ismember(cells.mouse, UI_mice))   = "UI";
    grp(ismember(cells.mouse, AVIR_mice)) = "Avirulent";
    P.group = categorical(grp);

    rng(42);
    P.Y = tsne(cells{:,1:cfg.nMetric}, 'Standardize', 1);

    P.types  = cfg.types;
    P.pretty = cfg.pretty;
end

function plot_by_celltype(P)
    hold on
    named = setdiff(categories(P.cellLabels), {'Other'}, 'stable');
    scatter(P.Y(P.cellLabels=="Other",1), P.Y(P.cellLabels=="Other",2), ...
        8, [0.85 0.85 0.85], 'filled', 'MarkerFaceAlpha', 0.3);
    cols = lines(numel(named));
    for i = 1:numel(named)
        idx = P.cellLabels==named{i};
        scatter(P.Y(idx,1), P.Y(idx,2), 12, cols(i,:), 'filled', 'MarkerFaceAlpha', 0.9);
    end
    
    leg = ['Other'; named(:)];
    for i = 1:numel(P.types)
        leg(strcmp(leg, P.types{i})) = P.pretty(i);
    end
    lgd = legend(leg, 'Location','bestoutside'); set(lgd,'Interpreter','tex','FontSize',10);
    xlabel('tSNE 1'); ylabel('tSNE 2'); set(gca,'FontName','Arial','FontSize',13); box on; hold off
end

function plot_by_group(P)
    hold on
    scatter(P.Y(P.group=="UI",1), P.Y(P.group=="UI",2), ...
        14, [0.2 0.6 1], 'filled', 'MarkerFaceAlpha', 0.7);
    scatter(P.Y(P.group=="Avirulent",1), P.Y(P.group=="Avirulent",2), ...
        14, [0.6 0.8 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
    legend({'Uninfected (UI)','ST1-75 colonized'}, 'Location','bestoutside','FontSize',10);
    xlabel('tSNE 1'); ylabel('tSNE 2'); set(gca,'FontName','Arial','FontSize',13); box on; hold off
end