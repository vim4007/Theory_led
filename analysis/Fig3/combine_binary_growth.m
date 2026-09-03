function combined = combine_binary_growth(parentFolder, outFile)
%
%   Usage:
%       M = combine_binary_growth('.');
%       M = combine_binary_growth('.', 'all_strains_growth.xlsx');

   dataDir = biologDataDir();                 % <repo>/data/Biolog
 
    if nargin < 1 || isempty(parentFolder) || parentFolder == ""
        parentFolder = fullfile(dataDir, 'Biolog_PM1_files');
    end
    if nargin < 2 || isempty(outFile) || outFile == ""
        outFile = fullfile(dataDir, 'Biolog_growth_matrix.xlsx');
    end
    parentFolder = char(parentFolder);
    outFile      = char(outFile);
 
    % Find strain-level binary growth files, dropping per-replicate ones.
    files = dir(fullfile(parentFolder, '**', '*_binary_growth.xlsx'));
    isRep = ~cellfun(@isempty, ...
        regexp({files.name}, '_binary_growth_\d+\.xlsx$', 'once'));
    files = files(~isRep);
    % Never re-ingest the combined output file itself.
    [~, outName, outExt] = fileparts(outFile);
    files = files(~strcmpi({files.name}, [outName outExt]));
    assert(~isempty(files), 'No *_binary_growth.xlsx files found under %s.', parentFolder);
 
    combined = table();
    for i = 1:numel(files)
        fpath  = fullfile(files(i).folder, files(i).name);
        strain = erase(string(files(i).name), "_binary_growth.xlsx");
        strainVar = matlab.lang.makeValidName(strain);
 
        t = readtable(fpath, 'VariableNamingRule', 'preserve');
        mets   = string(t{:, 1});
        growth = t{:, 2};
        gi = table(mets, growth, ...
                   'VariableNames', {'Metabolites', char(strainVar)});
 
        if isempty(combined)
            combined = gi;
        else
            combined = outerjoin(combined, gi, ...
                'Keys', 'Metabolites', 'MergeKeys', true);
        end
    end
 
    combined = sortrows(combined, 'Metabolites');
    writetable(combined, outFile);
    fprintf('Combined %d strains x %d metabolites -> %s\n', ...
        width(combined)-1, height(combined), outFile);
end
 