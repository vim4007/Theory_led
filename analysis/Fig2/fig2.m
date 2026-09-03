base_dir = fileparts(fileparts(fileparts(which('fig2'))));
vir_tbl = readtable(fullfile(base_dir, 'data', 'mouse', 'Scores', 'Virulence_screen_clean_table.csv'));
tbl     = readtable(fullfile(base_dir, 'data', 'mouse', 'Scores', 'ProtectionScreen_CDI_mouse.csv'));
tbl = [tbl; vir_tbl];
tbl.cdiffstrain = string(tbl.cdiffstrain);
tbl.tx          = string(tbl.tx);
tbl.experiment  = string(tbl.experiment);
tbl.exp_id      = tbl.experiment + "_" + tbl.cdiffstrain + "_" + tbl.mouse;

%%
no_vpi_st175_idx = ~contains(tbl.cdiffstrain, ".vpi") & ~contains(tbl.cdiffstrain, ".st1.75");
primary = tbl(no_vpi_st175_idx, :);                 

vpi_ending_idx = endsWith(tbl.cdiffstrain, ".vpi");
vpi_alone_idx  = strcmp(tbl.cdiffstrain, "vpi");
secondary = tbl(vpi_ending_idx | vpi_alone_idx, :); 

%% 
primary.experiment = categorical(primary.experiment);
primary.exp_id     = categorical(primary.exp_id);
primary.relweight(isnan(primary.relweight)) = 0;
primary.cdiffstrain = categorical(primary.cdiffstrain);
primary.cdiffstrain = reordercats(primary.cdiffstrain, ...
    ['ui'; setdiff(categories(primary.cdiffstrain), 'ui')]);
primary.relweight = -1 * primary.relweight; 

virulence_model = fitlme(primary, 'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
coeffs = fixedEffects(virulence_model);
names  = virulence_model.CoefficientNames;
ci     = coefCI(virulence_model);

virulence_score = table(names(:), coeffs(:), ci(:,1), ci(:,2), ...
    'VariableNames', {'cdiffstrain', 'Estimate', 'CI_Lower', 'CI_Upper'});
virulence_score = sortrows(virulence_score, 'Estimate', 'descend');
virulence_score.Strains = string(extractAfter(virulence_score.cdiffstrain, "cdiffstrain_"));

%%
secondary.experiment = categorical(secondary.experiment);
secondary.exp_id     = categorical(secondary.exp_id);
secondary.relweight(isnan(secondary.relweight)) = 0;
secondary.cdiffstrain = categorical(secondary.cdiffstrain);
secondary.cdiffstrain = reordercats(secondary.cdiffstrain, ...
    ['vpi'; setdiff(categories(secondary.cdiffstrain), 'vpi')]);

protection_model = fitlme(secondary, 'relweight ~ cdiffstrain + (1|day) + (1|exp_id)');
coeffs = fixedEffects(protection_model);
names  = protection_model.CoefficientNames;
ci     = coefCI(protection_model);

protection_score = table(names(:), coeffs(:), ci(:,1), ci(:,2), ...
    'VariableNames', {'cdiffstrain', 'Estimate', 'CI_Lower', 'CI_Upper'});
protection_score = sortrows(protection_score, 'Estimate', 'descend');
protection_score.Strains = extractAfter(protection_score.cdiffstrain, "cdiffstrain_");
protection_score.Strains = string(extractBefore(protection_score.Strains, ".vpi"));

%% 
combined_score = protection_score.Strains;
combined_score = array2table(combined_score);
combined_score(1,:) = []; 
combined_score = renamevars(combined_score, "combined_score", "Strains");

[lia, locb] = ismember(combined_score.Strains, virulence_score.Strains);
combined_score.Virulence_Estimate = NaN(height(combined_score), 1);
combined_score.Virulence_CI_Lower = NaN(height(combined_score), 1);
combined_score.Virulence_CI_Upper = NaN(height(combined_score), 1);
combined_score.Virulence_Estimate(lia) = virulence_score.Estimate(locb(lia));
combined_score.Virulence_CI_Lower(lia) = virulence_score.CI_Lower(locb(lia));
combined_score.Virulence_CI_Upper(lia) = virulence_score.CI_Upper(locb(lia));

[lia, locb] = ismember(combined_score.Strains, protection_score.Strains);
combined_score.Protection_Estimate = NaN(height(combined_score), 1);
combined_score.Protection_CI_Lower = NaN(height(combined_score), 1);
combined_score.Protection_CI_Upper = NaN(height(combined_score), 1);
combined_score.Protection_Estimate(lia) = protection_score.Estimate(locb(lia));
combined_score.Protection_CI_Lower(lia) = protection_score.CI_Lower(locb(lia));
combined_score.Protection_CI_Upper(lia) = protection_score.CI_Upper(locb(lia));

combined_score.Strains = upper(strrep(combined_score.Strains, '.', '-'));
combined_score(strcmp(combined_score.Strains, 'CD196'), :) = [];

arms      = ["vpi", "st1.75.vpi", "st1.68.vpi"];
armLabels = ["VPI", "ST1-75 + VPI", "ST1-68 + VPI"];
armColors = [0.85 0.15 0.15;   
             0.00 0.45 0.85;   
             0.00 0.60 0.30];  

D = tbl;
D.cdiffstrain = string(D.cdiffstrain);
D.exp_id      = string(D.exp_id);
D.relweight(isnan(D.relweight)) = 0;

%% 
figure('Color','w','Position',[100 100 520 400]); hold on
for a = 1:numel(arms)
    ad = D(D.cdiffstrain == arms(a), :);

    ids  = unique(ad.exp_id);
    keep = false(numel(ids),1);
    for i = 1:numel(ids)
        md = ad(ad.exp_id == ids(i), :);
        keep(i) = ~(any(md.death == 1) || any(md.relweight == 0));
    end
    surv = ad(ismember(ad.exp_id, ids(keep)), :);

    days = sort(unique(surv.day));
    m = NaN(numel(days),1); s = NaN(numel(days),1);
    for i = 1:numel(days)
        dm = surv.day == days(i);
        ex = surv.exp_id(dm); rw = surv.relweight(dm);
        ue = unique(ex);
        em = arrayfun(@(u) mean(rw(ex==u),'omitnan'), ue);
        em = em(~isnan(em));
        if ~isempty(em)
            m(i) = mean(em);
            if numel(em) > 1, s(i) = std(em); end
        end
    end
    v = ~isnan(m); x = days(v); y = m(v); e = s(v); e(isnan(e)) = 0;

    fill([x; flipud(x)], [y-e; flipud(y+e)], armColors(a,:), ...
        'FaceAlpha', 0.2, 'EdgeColor','none', 'HandleVisibility','off');
    plot(x, y, '-o', 'Color', armColors(a,:), 'MarkerFaceColor', armColors(a,:), ...
        'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', armLabels(a));
end
yline(100, 'k:', 'HandleVisibility','off');
xlabel('Day'); ylabel('Relative weight (%)');
legend('Location','southwest','Box','off');
xlim([0 7]); ylim([70 110]); box on; hold off
set(gca,'FontSize',25);
title('Co-colonization weight trajectories','FontSize',15);
%% 
figure('Color','w','Position',[100 100 520 400]); hold on
for a = 1:numel(arms)
    ad = D(D.cdiffstrain == arms(a), :);
    [G, ids] = findgroups(ad.exp_id);
    final_time  = splitapply(@max, ad.day,   G);
    final_death = splitapply(@max, ad.death, G);

    if sum(final_death) == 0
        x = [0; sort(unique(final_time)); max(final_time)];
        f = ones(size(x));
    else
        [f, x] = ecdf(final_time, 'Censoring', ~final_death, 'Function','survivor');
        if x(end) < max(final_time)
            x = [x; max(final_time)]; f = [f; f(end)];
        end
    end
    stairs(x, f, 'Color', armColors(a,:), 'LineWidth', 2.5, 'DisplayName', armLabels(a));
end
xlabel('Day'); ylabel('Survival probability');

legend('Location','southwest','Box','off');
ylim([0 1.05]); xlim([0 7]); box on; hold off
set(gca,'FontSize',25);
title('Co-colonization survival','FontSize',15);
%% 
cs = sortrows(combined_score, 'Protection_Estimate', 'descend');
err_lower = cs.Protection_Estimate - cs.Protection_CI_Lower;
err_upper = cs.Protection_CI_Upper - cs.Protection_Estimate;

figure('Color','w','Position',[100 100 1000 600]);
bar(cs.Protection_Estimate, 'FaceColor', [0.2 0.45 0.8]);
hold on;
errorbar(1:height(cs), cs.Protection_Estimate, err_lower, err_upper, ...
    'k.', 'LineWidth', 1.5);
hold off;
set(gca, 'XTick', 1:height(cs), 'XTickLabel', cs.Strains, 'XTickLabelRotation', 45);
xlabel('Strains'); ylabel('Protection score');
ylim([-10 25]); set(gca, 'FontSize', 18); grid on; box on;
set(gca,'FontSize',25);
title('Protection differs among ST1 strains','FontSize',15);
%%
cdt = combined_score;

hl = zeros(height(cdt),1);
hl(strcmp(cdt.Strains, 'ST1-75')) = 1;
hl(strcmp(cdt.Strains, 'ST1-68')) = 2;

blue  = [0 0.45 0.85];
green = [0 0.60 0.30];
gray  = [0.6 0.6 0.6];

figure('Color','w','Position',[100 100 800 600]); hold on

for i = 1:height(cdt)
    vl = cdt.Virulence_Estimate(i) - cdt.Virulence_CI_Lower(i);
    vu = cdt.Virulence_CI_Upper(i) - cdt.Virulence_Estimate(i);
    pl = cdt.Protection_Estimate(i) - cdt.Protection_CI_Lower(i);
    pu = cdt.Protection_CI_Upper(i) - cdt.Protection_Estimate(i);
    errorbar(cdt.Virulence_Estimate(i), cdt.Protection_Estimate(i), ...
        pl, pu, vl, vu, 'Color', [0.88 0.88 0.88], 'LineWidth', 1.5, ...
        'CapSize', 3, 'HandleVisibility','off');
end

mdl = fitlm(cdt, 'Protection_Estimate ~ Virulence_Estimate');
xr  = linspace(min(cdt.Virulence_Estimate), max(cdt.Virulence_Estimate), 100)';
yf  = predict(mdl, table(xr, 'VariableNames', {'Virulence_Estimate'}));
plot(xr, yf, 'k-', 'LineWidth', 3, 'HandleVisibility','off');

g = hl == 0;
scatter(cdt.Virulence_Estimate(g), cdt.Protection_Estimate(g), ...
    150, gray, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
scatter(cdt.Virulence_Estimate(hl==1), cdt.Protection_Estimate(hl==1), ...
    250, blue,  'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
scatter(cdt.Virulence_Estimate(hl==2), cdt.Protection_Estimate(hl==2), ...
    250, green, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

baseFont = 15;
for i = 1:height(cdt)
    if hl(i) == 0
        text(cdt.Virulence_Estimate(i)+0.2, cdt.Protection_Estimate(i)+0.2, ...
            cdt.Strains(i), 'FontSize', baseFont, 'Color', [0.35 0.35 0.35], ...
            'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    else
        c = blue; if hl(i)==2, c = green; end
        text(cdt.Virulence_Estimate(i)+0.25, cdt.Protection_Estimate(i)+0.25, ...
            cdt.Strains(i), 'FontSize', baseFont+10, 'FontWeight','bold', 'Color', c, ...
            'HorizontalAlignment','left', 'VerticalAlignment','bottom');
    end
end

[rho, pval] = corr(cdt.Virulence_Estimate, cdt.Protection_Estimate, ...
    'Type', 'Spearman', 'Rows', 'complete');
xl = xlim; yl = ylim;
text(xl(2), yl(2), sprintf('\\rho = %.2f, p = %.4f', rho, pval), ...
    'HorizontalAlignment','right', 'VerticalAlignment','top', 'FontSize', 16);

xlabel('Mono-colonization virulence score');
ylabel('Co-colonization protection score');

box on; hold off

[rP, pP] = corr(cdt.Virulence_Estimate, cdt.Protection_Estimate, 'Rows','complete');
fprintf('Spearman: rho=%.4f p=%.4g | Pearson: r=%.4f p=%.4g\n', rho, pval, rP, pP);
set(gca,'FontSize',25);
title('Low virulence alone does not ensure protection','FontSize',15);