
base_dir = fileparts(fileparts(fileparts(which('secondary_challnege'))));
tbl = readtable(fullfile(base_dir, 'data', 'mouse', 'Scores', 'weights.xlsx'));

%%
tbl.Date = datetime(tbl.Date);
first_date = min(tbl.Date);
tbl.Day = days(tbl.Date - first_date);

tbl.Group = cell(height(tbl), 1);
for i = 1:height(tbl)
    switch tbl.ExperimentalGroup{i}(1)
        case 'A', tbl.Group{i} = 'ST1-75';
        case 'B', tbl.Group{i} = 'ST1-68';
        case 'C', tbl.Group{i} = 'VPI';
    end
end

tbl.relweight = nan(height(tbl), 1);
unique_mice = unique(tbl.ExperimentalGroup);
for i = 1:length(unique_mice)
    mouse_idx = strcmp(tbl.ExperimentalGroup, unique_mice{i});
    mouse_data = tbl(mouse_idx, :);
    day0_weight = mouse_data.Weight_g_(mouse_data.Day == 0);
    tbl.relweight(mouse_idx) = (mouse_data.Weight_g_ / day0_weight) * 100;
end

tbl.death = double(isnan(tbl.relweight));
tbl.relweight(isnan(tbl.relweight)) = 0;

groups        = {'ST1-75', 'ST1-68', 'VPI'};
panel_titles  = {'ST1-75', 'ST1-68', 'VPI10463'};
group_colors  = {[0.85 0.20 0.20], [0.20 0.60 0.25], [0.15 0.35 0.75]};
weight_line   = [0 0.20 0.55]; 
weight_band   = [0.45 0.65 0.90]; 
challenge_days = [6 23];
fs = 20;

%%
figure('Position', [100 100 1500 450]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for g = 1:length(groups)
    gdata = tbl(strcmp(tbl.Group, groups{g}), :);
    mice  = unique(gdata.ExperimentalGroup);
    days_g = sort(unique(gdata.Day));

    mu = nan(length(days_g), 1);
    sd = nan(length(days_g), 1);

    for i = 1:length(days_g)
        vals = [];
        for m = 1:length(mice)
            md = sortrows(gdata(strcmp(gdata.ExperimentalGroup, mice{m}), :), 'Day');
            death_idx = find(md.death == 1, 1, 'first');
            if isempty(death_idx)
                death_day = Inf;
            else
                death_day = md.Day(death_idx);
            end
            day_idx = md.Day == days_g(i);
            if any(day_idx) && days_g(i) < death_day
                vals = [vals; md.relweight(day_idx)];
            end
        end
        if ~isempty(vals)
            mu(i) = mean(vals);
            if length(vals) > 1
                sd(i) = std(vals);
            else
                sd(i) = 0;
            end
        end
    end

    nexttile;
    hold on;

    for m = 1:length(mice)
        md = sortrows(gdata(strcmp(gdata.ExperimentalGroup, mice{m}), :), 'Day');
        death_idx = find(md.death == 1, 1, 'first');
        if ~isempty(death_idx) && death_idx > 1
            days_plot = md.Day(1:death_idx);
            w_plot = md.relweight(1:death_idx-1);
            w_plot = [w_plot; w_plot(end)];
            plot(days_plot, w_plot, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1);
            plot(md.Day(death_idx), w_plot(end), 'x', 'Color', [0.8 0 0], ...
                 'MarkerSize', 10, 'LineWidth', 2);
        end
    end

    valid = ~isnan(mu);
    x = days_g(valid);
    y = mu(valid);
    e = sd(valid);
    e(isnan(e)) = 0;
    fill([x; flipud(x)], [y-e; flipud(y+e)], weight_band, ...
         'FaceAlpha', 0.35, 'EdgeColor', 'none');
    plot(x, y, '-', 'Color', weight_line, 'LineWidth', 2);

    for d = challenge_days
        xline(d, 'k--', 'LineWidth', 1.2);
    end

    xlim([0 34]);
    ylim([80 115]);
    title(panel_titles{g}, 'FontSize', fs + 2);
    set(gca, 'FontSize', fs);
    grid on;
    box on;
    hold off;
end

xlabel(t, 'Day', 'FontSize', fs + 2);
ylabel(t, 'Relative weight (%)', 'FontSize', fs + 2);

%% 
figure('Position', [100 100 600 450]);
hold on;

for g = 1:length(groups)
    gdata = tbl(strcmp(tbl.Group, groups{g}), :);
    mice  = unique(gdata.ExperimentalGroup);
    days_g = sort(unique(gdata.Day));

    survival_prop = nan(length(days_g), 1);
    for i = 1:length(days_g)
        n_alive = 0;
        for m = 1:length(mice)
            md = sortrows(gdata(strcmp(gdata.ExperimentalGroup, mice{m}), :), 'Day');
            death_idx = find(md.death == 1, 1, 'first');
            if isempty(death_idx) || days_g(i) < md.Day(death_idx)
                n_alive = n_alive + 1;
            end
        end
        survival_prop(i) = n_alive / length(mice) * 100;
    end

    stairs(days_g, survival_prop, 'Color', group_colors{g}, 'LineWidth', 2, ...
           'DisplayName', panel_titles{g});
end

for d = challenge_days
    xline(d, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

xlabel('Day', 'FontSize', fs + 2);
ylabel('Survival (%)', 'FontSize', fs + 2);
xlim([0 35]);
ylim([0 105]);
legend('Location', 'southwest', 'FontSize', fs - 2);
grid on;
box on;
set(gca, 'FontSize', fs);
hold off;