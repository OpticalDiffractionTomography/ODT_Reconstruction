function p = phasePipeline(phaseIn, diskSE, dilSE, alreadyUnwrapped)
% phasePipeline  Wrapped phase -> flattened, aberration-compensated phase.
%
% Shared by field_Retrievals parallel pass (which passes wrapped phase) and
% its serial watchdog retry pass (which passes an already-unwrapped map, so
% unwrap2 is not called twice).
%
%   phaseIn           wrapped phase, or unwrapped if alreadyUnwrapped is true
%   diskSE, dilSE     strel objects for the cell-exclusion mask
%   alreadyUnwrapped  optional, default false

    if nargin < 4 || isempty(alreadyUnwrapped)
        alreadyUnwrapped = false;
    end

    if alreadyUnwrapped
        p = PhiShift(phaseIn);
    else
        p = PhiShift(unwrap2(phaseIn));
    end

    p = phaseCompensation(p, 1);

    % Cell-exclusion mask via top-hat filtering, so the degree-2 fit sees
    % background only.
    pp = p;
    tempThresh = p - imtophat(pp, diskSE);
    tempThresh = mean(tempThresh(:)) + 1.5;
    p2Mask = (pp > tempThresh(1));
    p2Mask = bwareaopen(p2Mask, 100);
    p2Mask = ~imdilate(p2Mask, dilSE);

    p = phaseCompensation(p, 2, p2Mask);
end
