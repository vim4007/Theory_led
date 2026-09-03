%run_biolog_folders();        % writes per-strain *_binary_growth.xlsx files
%combine_binary_growth();   % writes data/Biolog/Biolog_growth_matrix.xlsx

base_dir = fileparts(fileparts(fileparts(which('Fig3'))));
B   = readtable(fullfile(base_dir, 'data', 'Biolog', 'Biolog_growth_matrix.xlsx'), 'VariableNamingRule', 'preserve');
G   = readtable(fullfile(base_dir, 'data', 'Biolog', 'strain_groups.xlsx'),          'VariableNamingRule', 'preserve');
MoA = readtable(fullfile(base_dir, 'data', 'Biolog', 'moas.xlsx'),                    'VariableNamingRule', 'preserve');

%%
B = B(~strcmp(B.Metabolites,'Negative Control'),:);
G.Strains_norm = strrep(G.Strains,'-','_');
G_ST1 = G(startsWith(G.Strains_norm,'ST1'),:);
st1   = G_ST1.Strains_norm;
M     = table2array(B(:,st1));
prot  = G_ST1.Protection_Estimate;
n_met = height(B);  n_st1 = numel(st1);

MoA.category = regexprep(MoA.MoA,'^C-Source, ','');
[~, loc] = ismember(B.Metabolites, MoA.Chemical);
cats_per_row = MoA.category(loc);
cats = unique(cats_per_row(~cellfun(@isempty,cats_per_row)));

%% 
vpi_use = double(table2array(B(:,'VPI'))) > 0;    
[prot_sorted, ord] = sort(prot, 'descend');       
st1_sorted = st1(ord);
M_sorted   = M(:, ord);                            

hl_names = {'ST1_75','ST1_68','ST1_49'};
hl_col   = containers.Map({'ST1_75','ST1_68','ST1_49'}, ...
                          {[0 0.45 0.85],[0 0.60 0.30],[0.95 0.55 0.10]}); 

%% 
prev = sum(M,2);                                 
[~, ~, cat_idx] = unique(cats_per_row);            
cat_idx(cellfun(@isempty,cats_per_row)) = max(cat_idx)+1;  

sortkey = [~vpi_use, cat_idx, -prev];              
[~, col_ord] = sortrows(sortkey);

A_mat = [ vpi_use(col_ord)' ; M_sorted(col_ord,:)' ];   
row_labels = ['VPI10463'; strrep(st1_sorted,'ST1_','ST1-')];

figure('Color','w','Position',[100 100 900 620]);
imagesc(A_mat);
colormap([1 1 1; 0 0.30 0.75]); 
ax = gca;
set(ax,'FontSize',14);
lblFont = 2 * ax.FontSize;                    
set(ax,'YTick',1:numel(row_labels),'YTickLabel',row_labels);
set(ax,'XTick',[]);
xlabel('BIOLOG carbon sources','FontSize',lblFont);
ylabel('Strains sorted by protection','FontSize',lblFont);
title('Substrate-use matrix organized by VPI10463 use','FontSize',lblFont/2);
hold on;
n_vpi_used = sum(vpi_use);
xline(n_vpi_used + 0.5, 'k-', 'LineWidth', 1.5);
for k = 1:numel(hl_names)
    r = find(strcmp(st1_sorted, hl_names{k})) + 1;  
    if ~isempty(r)
        rectangle('Position',[0.5, r-0.5, n_met, 1], ...
            'EdgeColor',[0.95 0.55 0.10],'LineWidth',2);
    end
end


h_growth   = patch(NaN, NaN, [0 0.30 0.75], 'EdgeColor','k');
h_nogrowth = patch(NaN, NaN, [1 1 1],        'EdgeColor','k');
legend([h_growth, h_nogrowth], {'Growth','No growth'}, ...
    'Location','northeastoutside','Box','off','FontSize',lblFont);
hold off;
%% 
shared   = sum( M_sorted & vpi_use, 1)';
st1_priv = sum( M_sorted & ~vpi_use, 1)';
vpi_priv = sum(~M_sorted &  vpi_use, 1)';

figure('Color','w','Position',[100 100 950 500]);
bh = bar(1:n_st1, [shared, st1_priv, vpi_priv], 'stacked', 'EdgeColor','w');
bh(1).FaceColor = [0.6 0.6 0.6];
bh(2).FaceColor = [0 0.45 0.85];
bh(3).FaceColor = [0.85 0.15 0.15];
set(gca,'XTick',1:n_st1,'XTickLabel',strrep(st1_sorted,'ST1_','ST1-'), ...
    'XTickLabelRotation',45,'FontSize',12);
ylabel('Substrate count');
legend({'Shared','ST1 private','VPI private'},'Location','northeast','Box','off');
box off;
set(gca,'FontSize',20);
title('Pairwise resource classes','FontSize',15);
%%
breadth = sum(M,1)';

figure('Color','w','Position',[100 100 650 550]); hold on
scatter(breadth, prot, 150, [0.6 0.6 0.6], 'filled', 'MarkerEdgeColor','k');
p_fit = polyfit(breadth, prot, 1);
xr = linspace(min(breadth), max(breadth), 100);
plot(xr, polyval(p_fit, xr), 'k-', 'LineWidth', 2);
for k = 1:numel(hl_names)
    idx = strcmp(st1, hl_names{k});
    if any(idx)
        scatter(breadth(idx), prot(idx), 200, hl_col(hl_names{k}), 'filled', ...
            'MarkerEdgeColor','k','LineWidth',1.5);
        text(breadth(idx)-1.5, prot(idx)+0.7, strrep(hl_names{k},'ST1_','ST1-'), ...
            'FontWeight','bold','FontSize',20,'Color',hl_col(hl_names{k}));
    end
end
[rho_c, p_c_scatter] = corr(breadth, prot, 'Type','Spearman');
xl = xlim; yl = ylim;
text(xl(2), yl(1)+0.08*range(yl), sprintf('\\rho = %.2f, p = %.2f', rho_c, p_c_scatter), ...
    'HorizontalAlignment','right','FontSize',20);
xlabel('Total BIOLOG substrate-use breadth');
ylabel('Co-colonization Protection score');

set(gca,'FontSize',20); box on; hold off
title('Measured breadth vs protection','FontSize',15);
%% 
lam = linspace(0, 82, 501);
outcome_frac = zeros(n_st1, 4);   
for a = 1:n_st1
    s = M_sorted(:, a);
    nsh = sum( s &  vpi_use);
    nST1= sum( s & ~vpi_use);
    nVPI= sum(~s &  vpi_use);

    I_ST1 = nST1 - lam .* nVPI ./ (nsh + nVPI);
    I_VPI = nVPI - lam .* nST1 ./ (nsh + nST1);

    co = (I_ST1 >  0) & (I_VPI >  0);
    se = (I_ST1 >  0) & (I_VPI <= 0);
    ve = (I_ST1 <= 0) & (I_VPI >  0);
    ne = (I_ST1 <= 0) & (I_VPI <= 0);
    outcome_frac(a,:) = [mean(co), mean(se), mean(ve), mean(ne)];
end

figure('Color','w','Position',[100 100 950 500]);
bh = bar(1:n_st1, outcome_frac, 'stacked', 'EdgeColor','w');
bh(1).FaceColor = [0 0.60 0.30];
bh(2).FaceColor = [0 0.45 0.85];
bh(3).FaceColor = [0.85 0.15 0.15];
bh(4).FaceColor = [0.7 0.7 0.7];
set(gca,'XTick',1:n_st1,'XTickLabel',strrep(st1_sorted,'ST1_','ST1-'), ...
    'XTickLabelRotation',45,'FontSize',12);
ylabel('Fraction of \lambda grid'); ylim([0 1]);

legend({'Coexistence','ST1 excludes VPI','VPI excludes ST1','Neither Invades'}, ...
    'Location','northoutside','Orientation','horizontal','Box','off');
box off;
set(gca,'FontSize',20);
title('Model outcome classes','FontSize',15);
%% 
rho_e = nan(numel(cats),1); p_e = nan(numel(cats),1); tot_e = nan(numel(cats),1);
for j = 1:numel(cats)
    r = strcmp(cats_per_row, cats{j});
    tot_e(j) = sum(r);
    counts = sum(M(r,:),1)';
    if std(counts) > 0
        [rho_e(j), p_e(j)] = corr(counts, prot, 'Type','Spearman');
    end
end

q_e = nan(size(p_e));
valid = ~isnan(p_e);
pv = p_e(valid); [ps, si] = sort(pv); m = numel(ps);
qs = ps .* m ./ (1:m)';
qs = min(1, flipud(cummin(flipud(qs))));
qtmp = nan(size(pv)); qtmp(si) = qs;
q_e(valid) = qtmp;

E = table(cats, tot_e, rho_e, p_e, q_e, 'VariableNames', ...
    {'Category','TotalMets','rho','p','q'});
E = sortrows(E, 'rho', 'descend');

figure('Color','w','Position',[100 100 750 500]); hold on
bar(1:height(E), E.rho, 'FaceColor', [0.35 0.55 0.75], 'EdgeColor','k');
for i = 1:height(E)
    if E.q(i) < 0.05
        mk = '**';
    elseif E.p(i) < 0.05
        mk = '*';
    else
        mk = '';
    end
    if ~isempty(mk)
        yo = 0.03 * sign(E.rho(i)); va = 'bottom'; if E.rho(i)<0, va='top'; end
        text(i, E.rho(i)+yo, mk, 'FontSize',16,'FontWeight','bold', ...
            'HorizontalAlignment','center','VerticalAlignment',va);
    end
end
yline(0,'k');
set(gca,'XTick',1:height(E),'XTickLabel',E.Category,'XTickLabelRotation',45,'FontSize',13);
ylabel('Spearman \rho'); ylim([-0.2 0.8]);

box off; hold off

fprintf('\nPanel C: breadth vs protection  rho=%.2f p=%.3f\n', rho_c, p_c_scatter);
fprintf('Panel E categories (rho, p, q):\n');
for i = 1:height(E)
    fprintf('  %-18s rho=%+.2f  p=%.3f  q=%.3f\n', E.Category{i}, E.rho(i), E.p(i), E.q(i));
end
set(gca,'FontSize',20);
title('Chemical-class breadth vs protection','FontSize',15);