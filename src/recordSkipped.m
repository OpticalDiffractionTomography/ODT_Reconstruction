function recordSkipped(outDir, fileName, stage, reason)
% recordSkipped  Append a skipped/unreadable input to outDir/skipped_files.txt.
%   The orchestrator copies this file back with the results and reports it,
%   so a bad file is visible without having failed the whole chunk.
    fid = fopen(fullfile(outDir, 'skipped_files.txt'), 'a');
    if fid < 0, return; end
    fprintf(fid, '%s\t%s\t%s\n', fileName, stage, strrep(reason, newline, ' '));
    fclose(fid);
end
