function [Reconimg, ORytov] = ODTReconstruction(retAmplitude, retPhase, TP)

ORytov = gpuArray(single(zeros(TP.ZP2, TP.ZP2, TP.ZP3)));
Count  = single(zeros(TP.ZP2, TP.ZP2, TP.ZP3));

% Pre-compute the k-space grid once (identical for every frame)
kvec = TP.kres * (-floor(TP.ZP/2)+1 : floor(TP.ZP/2));
[ky, kx] = meshgrid(kvec, kvec);
kz = real(sqrt((TP.n_m*TP.k0)^2 - kx.^2 - ky.^2));
xr = TP.ZP * TP.res2 * TP.NA / TP.lambda;
NAmask  = ~mk_ellipse(xr, xr, TP.ZP, TP.ZP);  % logical, CPU
kzValid = (kz > 0) .* NAmask;                   % base validity mask

% Precompute k-volume bounds for fast index filtering
kxLo = TP.kres * (-floor(TP.ZP2/2)+1);
kyLo = TP.kres * (-floor(TP.ZP2/2)+1);
kzLo = TP.kres * (-floor(TP.ZP3/2)+1);
kxHi = TP.kres *  floor(TP.ZP2/2);
kyHi = TP.kres *  floor(TP.ZP2/2);
kzHi = TP.kres *  floor(TP.ZP3/2);

% Common padding amounts
padR = round((TP.ZP - TP.xx) / 2);
padC = round((TP.ZP - TP.yy) / 2);
res2sq = TP.res2^2;
fourPiKz = 1i * 4 * pi * kz;   % precomputed scalar field

for kk = TP.frameList
    FRytov = log(single(retAmplitude(:,:,kk))) + 1i*single(retPhase(:,:,kk));
    FRytov = gpuArray(padarray(FRytov, [padR padC], 'symmetric'));

    UsRytov = fftshift(fft2(FRytov)) * res2sq;
    UsRytov = circshift(UsRytov, [round(TP.f_dx2(kk)) round(TP.f_dy2(kk))]);
    UsRytov = UsRytov .* NAmask;

    Kx = kx - TP.k0_x(kk);
    Ky = ky - TP.k0_y(kk);
    Kz = kz - TP.k0_z(kk);
    Uprime = fourPiKz .* UsRytov;

    xind = find(kzValid ...
        .* (Kx > kxLo) .* (Kx < kxHi) ...
        .* (Ky > kyLo) .* (Ky < kyHi) ...
        .* (Kz > kzLo) .* (Kz < kzHi));

    Uprime = Uprime(xind);
    KindX  = round(Kx(xind) / TP.kres + TP.ZP2/2);
    KindY  = round(Ky(xind) / TP.kres + TP.ZP2/2);
    KindZ  = round(Kz(xind) / TP.kres + TP.ZP3/2);
    Klinear = (KindZ-1)*TP.ZP2^2 + (KindY-1)*TP.ZP2 + KindX;

    [~, idx] = unique(Klinear, 'first');
    idx      = sort(idx);
    Klinear  = Klinear(idx);
    Uprime   = Uprime(idx);

    ORytov(Klinear) = ORytov(Klinear) + Uprime;
    Count(Klinear)  = Count(Klinear) + 1;
end

ORytov = gather(ORytov);
ORytov(Count > 0) = ORytov(Count > 0) ./ Count(Count > 0);
ORytov = gpuArray(ORytov);

Reconimg = ifftn(ifftshift(ORytov)) / (TP.res3^2 * TP.res4);
Reconimg = TP.n_m * sqrt(1 - Reconimg .* (TP.lambda / (TP.n_m*2*pi))^2);
Reconimg = gather(Reconimg);
return;