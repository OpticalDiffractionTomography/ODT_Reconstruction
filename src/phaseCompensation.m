function [goodp2, coefficients, compensationMap] = phaseCompensation(varargin)
p2 = varargin{1};
n  = varargin{2};

if length(varargin) == 3
    mask = varargin{3};
else
    mask = ones(size(p2));
end

[imY, imX] = size(p2);
[XX, YY]   = meshgrid(1:imX, 1:imY);

% Build polynomial basis columns in one vectorized step (no loop)
% Powers 1..n for X, then 1..n for Y, then constant — matches original layout
pwrs = (1:n);
XXv  = XX(:);  YYv = YY(:);
Xbasis = XXv .^ pwrs;   % (imX*imY) x n
Ybasis = YYv .^ pwrs;   % (imX*imY) x n
AA_full = [Xbasis, Ybasis, ones(imX*imY, 1)];  % (imX*imY) x (2n+1)

% Apply mask: keep only unmasked pixels for the least-squares fit
p2vec = p2(:) .* mask(:);
keep  = p2vec ~= 0;
AA_masked  = AA_full(keep, :);
p2_masked  = p2vec(keep);

coefficients = (AA_masked' * AA_masked) \ (AA_masked' * p2_masked);

% Reconstruct compensation map over the full image (vectorized)
compensationMap = reshape(AA_full * coefficients, imY, imX);

goodp2 = p2 - compensationMap;