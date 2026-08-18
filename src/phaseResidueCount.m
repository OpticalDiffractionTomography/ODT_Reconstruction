function nRes = phaseResidueCount(wrappedPhase)
% phaseResidueCount  Number of phase residues in a wrapped phase map.
%
% A residue is a 2x2 loop whose wrapped phase differences sum to a non-zero
% multiple of 2*pi — i.e. a point where the phase field is not integrable.
% Goldstein branch-cut unwrapping must pair these up, and its cost grows
% steeply with their number, which is what makes strongly scattering samples
% hang. Counting them first is cheap (a few vectorized diffs) and lets
% field_Retrieval route risky frames to a watchdogged retry instead of
% gambling the whole job on them.

    w = @(d) mod(d + pi, 2*pi) - pi;   % wrap differences into (-pi, pi]

    % Closed loop around each 2x2 neighbourhood:
    %   (i,j) -> (i,j+1) -> (i+1,j+1) -> (i+1,j) -> (i,j)
    a = w(wrappedPhase(1:end-1, 2:end  ) - wrappedPhase(1:end-1, 1:end-1));
    b = w(wrappedPhase(2:end,   2:end  ) - wrappedPhase(1:end-1, 2:end  ));
    c = w(wrappedPhase(2:end,   1:end-1) - wrappedPhase(2:end,   2:end  ));
    d = w(wrappedPhase(1:end-1, 1:end-1) - wrappedPhase(2:end,   1:end-1));

    loopSum = a + b + c + d;

    % Non-zero loop sums are +/-2*pi (or multiples); tolerance guards rounding.
    nRes = sum(abs(loopSum(:)) > pi);
end
