function p = unwrapWithTimeout(wrappedPhase, timeoutSec)
% unwrapWithTimeout  Run unwrap2 with a hard wall-clock limit.
%
% unwrap2 (Goldstein branch-cut) can spin effectively forever on
% residue-heavy frames from strongly scattering samples. A try/catch cannot
% interrupt that — it never throws — so the call is pushed onto a parallel
% pool and abandoned if it overruns.
%
% Returns the unwrapped phase, or throws 'ODT:unwrapTimeout' so the caller
% can NaN-fill the frame and move on.
%
% Note: this must NOT be called from inside a parfor body (nested pool use is
% not allowed). field_Retrieval calls it from a serial retry pass.

    pool = gcp('nocreate');
    if isempty(pool)
        % No pool available: run inline, unprotected. Better than failing
        % outright, but the hang risk is back.
        p = unwrap2(wrappedPhase);
        return
    end

    f = parfeval(pool, @unwrap2, 1, wrappedPhase);
    ok = wait(f, 'finished', timeoutSec);

    if ~ok
        cancel(f);
        error('ODT:unwrapTimeout', ...
            'unwrap2 exceeded %g s and was abandoned', timeoutSec);
    end

    if ~isempty(f.Error)
        err = f.Error;
        cancel(f);
        rethrow(err);
    end

    p = fetchOutputs(f);
end
