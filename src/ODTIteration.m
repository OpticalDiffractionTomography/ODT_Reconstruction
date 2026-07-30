function Reconimg = ODTIteration(Reconimg, ORytov, TP, IterNum)

ORytov_index = find(abs(ORytov) > 0);
ORytov       = ORytov(ORytov_index);

% Precompute loop-invariant scalars on GPU
normFact  = 1 / (TP.res3^2 * TP.res4);
rytovScale = (TP.lambda / (TP.n_m * 2*pi))^2;  % used in sqrt step
chiScale   = -(2*pi*TP.n_m / TP.lambda)^2;     % forward scattering potential
nm_inv2    = 1 / TP.n_m^2;

for mm = 1:IterNum
    % Physical constraint: real(n) >= n_m
    id = real(Reconimg) < TP.n_m;
    Reconimg(id) = TP.n_m - 1i*imag(Reconimg(id));

    % Convert RI to scattering potential chi
    Reconimg = chiScale .* (Reconimg.^2 * nm_inv2 - 1);

    % Project onto measured Fourier voxels
    ORytov_new = fftshift(fftn(Reconimg)) / normFact;
    ORytov_new(ORytov_index) = ORytov;

    % Back to RI
    Reconimg = ifftn(ifftshift(ORytov_new)) * normFact;
    Reconimg = TP.n_m * sqrt(1 - Reconimg * rytovScale);
end

Reconimg = gather(fftshift(Reconimg, 3));
return