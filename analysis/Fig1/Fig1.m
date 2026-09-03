base_dir = fileparts(fileparts(fileparts(which('Fig1'))));
B   = readtable(fullfile(base_dir, 'data', 'Biolog', 'Biolog_growth_matrix.xlsx'), 'VariableNamingRule', 'preserve');
G   = readtable(fullfile(base_dir, 'data', 'Biolog', 'strain_groups.xlsx'), 'VariableNamingRule', 'preserve');
%%
vpi = B.VPI;

counts = @(s) deal( sum(s & vpi), ...      % n_shared
                    sum(s & ~vpi), ...     % n_ST1  (ST1-private)
                    sum(~s & vpi) );       % n_VPI  (VPI-private)

[nsh_75, nST1_75, nVPI_75] = counts(B.ST1_75);
[nsh_68, nST1_68, nVPI_68] = counts(B.ST1_68);

fprintf('ST1-75:  n_shared=%d  n_ST1=%d  n_VPI=%d  (shared+VPI=%d)\n', ...
        nsh_75, nST1_75, nVPI_75, nsh_75+nVPI_75);
fprintf('ST1-68:  n_shared=%d  n_ST1=%d  n_VPI=%d  (shared+VPI=%d)\n', ...
        nsh_68, nST1_68, nVPI_68, nsh_68+nVPI_68);

I_ST1 = @(nST1,nVPI,nsh,lam) nST1 - lam .* nVPI ./ (nsh + nVPI);
I_VPI = @(nVPI,nST1,nsh,lam) nVPI - lam .* nST1 ./ (nsh + nST1);

lam = linspace(0, 90, 200);

blue = [0 0.45 0.85];
red  = [0.85 0.15 0.15];

%% 
figure('Color','w','Position',[100 100 420 380]); hold on
L = 10;
fill([0 L L 0],[0 0 L L], [0.80 0.92 0.80],'EdgeColor','none'); % coexistence  (both +)
fill([-L 0 0 -L],[0 0 L L],[0.96 0.80 0.80],'EdgeColor','none'); % VPI excludes ST1
fill([0 L L 0],[-L -L 0 0],[0.80 0.86 0.96],'EdgeColor','none'); % ST1 excludes VPI
fill([-L 0 0 -L],[-L -L 0 0],[0.88 0.88 0.88],'EdgeColor','none');% neither
plot([-L L],[0 0],'k:','LineWidth',1);
plot([0 0],[-L L],'k:','LineWidth',1);
text(-L/2, L/2,  'VPI excludes ST1','HorizontalAlignment','center','FontSize',20);
text( L/2, L/2,  'Coexistence',      'HorizontalAlignment','center','FontSize',20);
text(-L/2,-L/2,  'Neither Invades',  'HorizontalAlignment','center','FontSize',20);
text( L/2,-L/2,  'ST1 excludes VPI', 'HorizontalAlignment','center','FontSize',20);
axis([-L L -L L]); axis square; box on
xlabel('ST1 invasion score  I_{ST1}');
ylabel('VPI invasion score  I_{VPI}');
set(gca,'FontSize',20);
title('Outcome rule from invasion scores');
hold off

%%
figure('Color','w','Position',[100 100 460 360]); hold on
plot(lam, I_ST1(nST1_75,nVPI_75,nsh_75,lam),'Color',blue,'LineWidth',4);
plot(lam, I_VPI(nVPI_75,nST1_75,nsh_75,lam),'Color',red ,'LineWidth',4);
yline(0,'k:');
xlabel('Resource limitation threshold  \lambda');
ylabel('Invasion score');
title('ST1-75 / VPI10463-like counts');
legend({'ST1 invades VPI','VPI invades ST1'},'Location','southwest','Box','off');
set(gca,'FontSize',20);
xlim([0 90]); ylim([-70 70]); box on; hold off

%% 
figure('Color','w','Position',[100 100 460 360]); hold on
plot(lam, I_ST1(nST1_68,nVPI_68,nsh_68,lam),'Color',[0 0.6 0.3],'LineWidth',4);
plot(lam, I_VPI(nVPI_68,nST1_68,nsh_68,lam),'Color',red ,'LineWidth',4);
yline(0,'k:');
xlabel('Resource limitation threshold  \lambda');
ylabel('Invasion score');
title('ST1-68 / VPI10463-like counts');
legend({'ST1 invades VPI','VPI invades ST1'},'Location','southwest','Box','off');
set(gca,'FontSize',20);
xlim([0 90]); ylim([-70 70]); box on; hold off