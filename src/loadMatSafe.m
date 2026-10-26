function [S, errMsg] = loadMatSafe(matPath, requiredVars)
% loadMatSafe  Load a MAT-file without letting a corrupt file abort the job.
%   [S, errMsg] = loadMatSafe(matPath, requiredVars) returns the loaded struct S
%   and errMsg = '' on success. If the file cannot be read, or lacks any of the
%   variables in requiredVars (cellstr), S = [] and errMsg describes why.
    S = [];
    errMsg = '';
    try
        S = load(matPath);
    catch ME
        S = [];
        errMsg = sprintf('cannot read MAT-file: %s', ME.message);
        return
    end
    missing = requiredVars(~isfield(S, requiredVars));
    if ~isempty(missing)
        S = [];
        errMsg = sprintf('missing variable(s): %s', strjoin(missing, ', '));
    end
end
