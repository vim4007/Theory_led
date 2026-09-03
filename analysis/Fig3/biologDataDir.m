function dataDir = biologDataDir()
% BIOLOGDATADIR  Return the absolute path to <repo>/data/Biolog.

    thisFile = mfilename('fullpath');            % .../code/.../biologDataDir
    repoRoot = fileparts(fileparts(fileparts(thisFile)));
    dataDir  = fullfile(repoRoot, 'data', 'Biolog');

    if ~isfolder(dataDir)
        error(['Could not locate data/Biolog relative to %s.\n' ...
               'Expected: %s\nCheck that this file sits three levels ' ...
               'above the repo root, or edit biologDataDir.m.'], ...
               thisFile, dataDir);
    end
end