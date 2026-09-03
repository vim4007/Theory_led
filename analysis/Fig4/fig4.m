
base_dir = fileparts(fileparts(fileparts(which('fig4'))));
B = readtable(fullfile(base_dir, 'data', 'Biolog', 'Biolog_growth_matrix.xlsx'), ...
    'VariableNamingRule', 'preserve');
B = B(~strcmp(B.Metabolites,'Negative Control'),:);

weights = readtable(fullfile(base_dir, 'data', 'qPCR', 'weights.xlsx'));
qpcr    = readtable(fullfile(base_dir, 'data', 'qPCR', 'qPCR section.xlsx'));
vpi_use = double(table2array(B(:,'VPI'))) > 0;

countsFor = @(name) deal( ...
    sum( table2array(B(:,name))>0 &  vpi_use), ...   % n_shared
    sum( table2array(B(:,name))>0 & ~vpi_use), ...   % n_ST1
    sum(~(table2array(B(:,name))>0) &  vpi_use));     % n_VPI

[nsh75, nST1_75, nVPI_75] = countsFor('ST1_75');
[nsh68, nST1_68, nVPI_68] = countsFor('ST1_68');

fprintf('ST1-75: shared=%d ST1-priv=%d VPI-priv=%d\n', nsh75, nST1_75, nVPI_75);
fprintf('ST1-68: shared=%d ST1-priv=%d VPI-priv=%d\n', nsh68, nST1_68, nVPI_68);

blue  = [0 0.45 0.85];
green = [0 0.60 0.30];
red   = [0.85 0.15 0.15];

%% 
lam = linspace(0, 82, 501);
I_ST1 = @(nST1,nVPI,nsh,l) nST1 - l .* nVPI ./ (nsh + nVPI);

figure('Color','w','Position',[100 100 560 440]); hold on
plot(lam, I_ST1(nST1_75,nVPI_75,nsh75,lam), 'Color', blue,  'LineWidth', 4.5, 'DisplayName','ST1-75');
plot(lam, I_ST1(nST1_68,nVPI_68,nsh68,lam), 'Color', green, 'LineWidth', 4.5, 'DisplayName','ST1-68');
yline(0,'k:','HandleVisibility','off');
xlabel('Loss/supply threshold  \lambda');
ylabel('ST1 invasion score');

legend('Location','southwest','Box','off');
xlim([0 90]); ylim([-40 80]); set(gca,'FontSize',20); box on; hold off
title('Predicted invasion robustness','FontSize',15);
%% 
e      = 0.5;
d      = 0.05;
D      = 0.2;
delta0 = 0.08;
u      = 0.25;

n_sh  = nsh75;
n_st1 = nST1_75;
n_vpi = nVPI_75;

odef = @(t,y) [ ...
    y(1) * ( e*u*y(3) + e*u*y(4) - d - D ); ...
    y(2) * ( e*u*y(3) + e*u*y(5) - d - D ); ...
    n_sh  * delta0 - D*y(3) - u*y(3)*(y(1)+y(2)); ...
    n_st1 * delta0 - D*y(4) - u*y(4)* y(1); ...
    n_vpi * delta0 - D*y(5) - u*y(5)* y(2) ];

N_ST1_0 = 1;  N_VPI_0 = 5;
y0 = [N_ST1_0; N_VPI_0; ...
      n_sh*delta0/D; n_st1*delta0/D; n_vpi*delta0/D];

tspan = [0 25];
[t, Y] = ode45(odef, tspan, y0);

figure('Color','w','Position',[100 100 620 440]); hold on
plot(t, Y(:,1), 'Color', blue, 'LineWidth', 4.5, 'DisplayName','ST1-75');
plot(t, Y(:,2), 'Color', red,  'LineWidth', 4.5, 'DisplayName','VPI');
xlabel('Time (model units)');
ylabel('Density');

legend('Location','east','Box','off');
xlim([0 25]); set(gca,'FontSize',20); box on; hold off

lambda_val = (d+D)*D / (e*u*delta0);
fprintf('Implied lambda = %.2f\n', lambda_val);
title('Model initialized at 1:5 ST1-75:VPI','FontSize',15);
%%
weights = renamevars(weights, "Var1", "Days");
blue = [0.20 0.55 0.80];
red  = [0.85 0.33 0.10];
gray = [0.4 0.4 0.4];

%% 
days_w   = weights.Days;
wt_st175 = weights.ST1_75;
wt_vpi   = weights.VPI;
wt_mix   = [weights.mix, weights.mix_1, weights.mix_2];

mix_mean = mean(wt_mix, 2);
mix_std  = std(wt_mix, 0, 2);

figure('Color','w','Position',[100 100 600 450]); hold on
fill([days_w; flipud(days_w)], [mix_mean - mix_std; flipud(mix_mean + mix_std)], ...
     [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor','none', 'HandleVisibility','off');

plot(days_w, wt_st175, '-o', 'Color', blue, 'LineWidth', 4, 'MarkerSize', 7, ...
    'MarkerFaceColor', blue, 'DisplayName','ST1-75');
plot(days_w, wt_vpi, '-o', 'Color', red, 'LineWidth', 4, 'MarkerSize', 7, ...
    'MarkerFaceColor', red, 'DisplayName','VPI');
plot(days_w, mix_mean, '-o', 'Color', gray, 'LineWidth', 4, 'MarkerSize', 7, ...
    'MarkerFaceColor', gray, 'DisplayName','ST1-75 + VPI');

yline(100,'k--','LineWidth',1,'HandleVisibility','off');
xlabel('Days post-infection');
ylabel('Relative weight (%)');

legend('Location','southwest','Box','off','FontSize',14);
xlim([0 3.2]); ylim([60 115]);
set(gca,'XTick',[0 1 2 3],'FontSize',20); box off; hold off
title('Weight after 1:5 ST1-75:VPI inoculation','FontSize',15);
%%
days_q = [0.5, 1, 3];
mean_co_st175 = nan(3,1);
mean_co_vpi   = nan(3,1);
for d = 1:3
    day_data = qpcr(qpcr.Day == days_q(d), :);
    co_idx = ismember(day_data.Mouse, [3,4,5]);
    mean_co_st175(d) = nanmean(day_data.ST1_75(co_idx));
    mean_co_vpi(d)   = nanmean(day_data.VPI(co_idx));
end

figure('Color','w','Position',[100 100 560 450]);
b = bar(1:3, [mean_co_vpi, mean_co_st175], 'stacked', 'EdgeColor','k','LineWidth',1.2);
b(1).FaceColor = [1 0 0];      
b(2).FaceColor = [0 0 0.8];   
set(gca,'XTick',1:3,'XTickLabel',{'0.5','1','3'},'FontSize',16);
xlabel('Day');
ylabel('Mean qPCR fraction (%)');

ylim([0 100]); xlim([0.4 3.6]);
legend({'VPI','ST1-75'},'Location','northwest','Box','on','FontSize',20);
box on; set(gca,'LineWidth',1.2);

fprintf('Co-infection ST1-75 fraction: day0.5=%.0f%% day1=%.0f%% day3=%.0f%%\n', ...
    mean_co_st175(1), mean_co_st175(2), mean_co_st175(3));
set(gca,'FontSize',20);
title('ST1-75 increases during co-infection','FontSize',15);