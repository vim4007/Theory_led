function run_biolog_folders(parentFolder, biologNamesFile, opts)
%   
%   Example (run from the parent folder that holds analyze_biolog_plate.m):
%       run_biolog_folders('.', 'Biolog_names_sorted.xlsx');
%       run_biolog_folders('.', 'Biolog_names_sorted.xlsx', 'Plot', true);

      arguments
        parentFolder      (1,1) string = ""
        biologNamesFile   (1,1) string = ""
        opts.Plot         (1,1) logical = false
        opts.TimeCutoff   (1,1) double  = 24
        opts.Alpha        (1,1) double  = 0.05
        opts.NumTests     (1,1) double  = 95
        opts.Pattern      (1,1) string  = "*_PM1*.xlsx"
    end
 
    dataDir = biologDataDir();                 % <repo>/data/Biolog
    if parentFolder == ""
        parentFolder = fullfile(dataDir, 'Biolog_PM1_files');
    end
    if biologNamesFile == ""
        biologNamesFile = fullfile(parentFolder, 'Biolog_names_sorted.xlsx');
    end
    assert(isfolder(parentFolder), 'Plate folder not found: %s', parentFolder);
    assert(isfile(biologNamesFile), 'Biolog names file not found: %s', biologNamesFile);
    biologNamesFile = string(biologNamesFile);
 
    plateOpts = {'Plot', opts.Plot, ...
                 'TimeCutoff', opts.TimeCutoff, ...
                 'Alpha', opts.Alpha, 'NumTests', opts.NumTests};
 
    d = dir(parentFolder);
    d = d([d.isdir]);
    d = d(~ismember({d.name}, {'.', '..'}));
    assert(~isempty(d), 'No subfolders found in %s.', parentFolder);
 
    for s = 1:numel(d)
        strainDir  = fullfile(d(s).folder, d(s).name);
        strainName = string(d(s).name);
 
        files = dir(fullfile(strainDir, opts.Pattern));
        if isempty(files)
            warning('%s: no files matching %s, skipping.', strainName, opts.Pattern);
            continue;
        end
        fpaths = fullfile(strainDir, string({files.name}'));
        fpaths = sortNat(fpaths);
 
        fprintf('\n=== %s (%d file(s)) ===\n', strainName, numel(fpaths));
 
        if numel(fpaths) == 1
            g = analyze_biolog_plate(fpaths(1), biologNamesFile, plateOpts{:});
            outFile = fullfile(strainDir, strainName + "_binary_growth.xlsx");
            writetable(g, outFile);
            fprintf('Wrote %s\n', outFile);
        else
            combined = table();
            for i = 1:numel(fpaths)
                g = analyze_biolog_plate(fpaths(i), biologNamesFile, plateOpts{:});
 
                [~, repBase] = fileparts(fpaths(i));
                repTok = regexp(repBase, '_PM\d+_(\d+)$', 'tokens', 'once');
                if isempty(repTok); repNum = i; else; repNum = str2double(repTok{1}); end
 
                repFile = fullfile(strainDir, ...
                    strainName + "_binary_growth_" + repNum + ".xlsx");
                writetable(g, repFile);
 
                gi = g(:, {'Metabolites'});
                gi.("rep" + i) = g.Growth;
                if isempty(combined)
                    combined = gi;
                else
                    combined = outerjoin(combined, gi, ...
                        'Keys', 'Metabolites', 'MergeKeys', true);
                end
            end
 
            repData = combined{:, 2:end};
            repData(isnan(repData)) = 0;
            final = table();
            final.Metabolites = combined.Metabolites;
            final.Growth      = double(any(repData == 1, 2));
 
            outFile = fullfile(strainDir, strainName + "_binary_growth.xlsx");
            writetable(final, outFile);
            fprintf('%s: OR of %d replicates -> %d / %d support growth.\n', ...
                strainName, numel(fpaths), sum(final.Growth), height(final));
            fprintf('Wrote %s\n', outFile);
        end
    end
end
 
function s = sortNat(s)
    n = numel(s);
    key = zeros(n, 1);
    for i = 1:n
        [~, b] = fileparts(s(i));
        d = regexp(b, '_(\d+)$', 'tokens', 'once');
        if isempty(d); key(i) = 0; else; key(i) = str2double(d{1}); end
    end
    [~, ord] = sort(key);
    s = s(ord);
end
 
