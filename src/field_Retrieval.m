%% Field Retrieval
% spath is provided by main.sh at invocation
set(0, 'DefaultFigureVisible', 'off');  % headless: no display on compute nodes

outDir = fullfile(spath, 'field_retrieval');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logFile = fullfile(outDir, 'field_retrieval.log');
fid = fopen(logFile, 'w');
logfn = @(msg) fprintf(fid, '[%s] %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), msg);
logfn('=== field_Retrieval started ===');
logfn(sprintf('spath: %s', spath));

%% Start parallel pool (workers = SLURM CPUs minus 1 for main thread)
nWorkers = max(1, str2double(getenv('SLURM_CPUS_PER_TASK')) - 1);
if isnan(nWorkers) || nWorkers < 1, nWorkers = 4; end
if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end
logfn(sprintf('Parallel pool: %d workers', nWorkers));

%% Background loading
logfn('Scanning for background files...');
bglist = dir(fullfile(spath, 'bg*_Tomog.mat'));
logfn(sprintf('Found %d background file(s). Loading: %s', length(bglist), bglist(1).name));
f_dx = [];
f_dy = [];
load(fullfile(spath, bglist(1).name));
img = double(squeeze(tomogMap(:,:,round(size(tomogMap,3)/3)-1)));
img = squeeze(img);
ii = length(img);
ZP = ii;
r = round(ZP*res*NA/lambda)+20;
yr = r;
logfn(sprintf('Image size: %d px | ZP=%d | r=%d | NA=%.3f | lambda=%.4f um | res=%.4f um', ii, ZP, r, NA, lambda, res));

%% Carrier frequency detection from background frames
logfn(sprintf('Detecting carrier frequency from %d background frames...', size(tomogMap,3)));
nBg = size(tomogMap,3);
f_dx_bg = zeros(nBg,1);
f_dy_bg = zeros(nBg,1);
bgStack = double(tomogMap);  % load once; avoids repeated squeeze+cast inside loop
parfor bgNum = 1:nBg
    bgImg = bgStack(:,:,bgNum);
    Ftmp = fftshift(fft2(bgImg))/(ii^2);
    [fx,fy] = find(Ftmp==max(max(Ftmp(:, round(ii*0.01):round(ii*0.49)))));
    f_dx_bg(bgNum) = fx(1);
    f_dy_bg(bgNum) = fy(1);
end
f_dx = f_dx_bg;
f_dy = f_dy_bg;
mi = mean(f_dx);
mj = mean(f_dy);
mi = round(mi-ii/2-1); mj = round(mj-ii/2-1);
logfn(sprintf('Carrier frequency offset: mi=%d px, mj=%d px', mi, mj));

%% Build demodulation mask and background field stack
logfn('Building demodulation mask (mk_ellipse)...');
c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
c3mask = circshift(c1mask,[mi mj]);
logfn(sprintf('Demodulating %d background frames...', size(tomogMap,3)));
nBg    = size(tomogMap,3);
Fbg    = zeros(round(2*r),round(2*r),nBg,'single');
bgStack = double(tomogMap);  % cast once
mi_r = round(mi); mj_r = round(mj);
rowRange = ii/2-r+1:ii/2+r;
parfor bgNum = 1:nBg
    bgImg = bgStack(:,:,bgNum);
    Ftmp = fftshift(fft2(bgImg))/(ii^2);
    Ftmp = Ftmp.*c3mask;
    Ftmp = circshift(Ftmp,-[mi_r mj_r]);
    Ftmp = Ftmp(rowRange, rowRange);
    sz   = length(Ftmp);
    Fbg(:,:,bgNum) = single(ifft2(ifftshift(Ftmp))*(sz^2));
end
logfn('Background field stack ready.');

%% Per-sample field retrieval
sampleList = dir(fullfile(spath, 'sample*_Tomog.mat'));
logfn(sprintf('Found %d sample file(s).', length(sampleList)));

for sampleNum = 1:length(sampleList)
    sName = sampleList(sampleNum).name;
    logfn(sprintf('--- Processing sample %d/%d: %s ---', sampleNum, length(sampleList), sName));

    logfn(sprintf('  Loading %s...', sName));
    load(fullfile(spath, sName));
    nFrames = size(tomogMap,3);
    logfn(sprintf('  Loaded: %d frames, size %dx%d', nFrames, size(tomogMap,1), size(tomogMap,2)));

    pSize = round(2*r)-4;
    retPhase     = zeros(pSize, pSize, nFrames, 'single');
    retAmplitude = zeros(pSize, pSize, nFrames, 'single');
    f_dx_s = zeros(nFrames,1);
    f_dy_s = zeros(nFrames,1);

    % Pre-cast the whole stack once to avoid repeated double() inside parfor
    sampleStack = double(tomogMap);

    mi_r = round(mi); mj_r = round(mj);
    rowRange = ii/2-r+1:ii/2+r;
    ii_sq  = ii^2;
    colLo  = round(ii*0.01);
    colHi  = round(ii*0.49);
    diskSE = strel('disk', 150);  % created once; safe to broadcast (value object)
    dilSE  = strel('disk', 5);

    logfn(sprintf('  Launching parfor over %d frames (%d workers)...', nFrames, nWorkers));
    parfor iter = 1:nFrames
        frImg = sampleStack(:,:,iter);
        Fimg  = fftshift(fft2(frImg)) / ii_sq;

        % carrier frequency detection
        [f_x,f_y] = find(Fimg == max(max(Fimg(:, colLo:colHi))));
        f_dx_s(iter) = f_x(1);
        f_dy_s(iter) = f_y(1);

        % demodulation
        Fimg = Fimg .* c3mask;
        Fimg = circshift(Fimg, -[mi_r mj_r]);
        Fimg = Fimg(rowRange, rowRange);
        sz   = length(Fimg);
        Fimg = ifft2(ifftshift(Fimg)) * (sz^2);
        Fimg = Fimg ./ Fbg(:,:,iter);

        % residual tilt correction via peak shift
        FFimg = fftshift(fft2(Fimg));
        [tX,tY] = find(abs(FFimg) == max(max(abs(FFimg))));
        tX = tX(1) - size(FFimg,1)/2;
        tY = tY(1) - size(FFimg,2)/2;
        Fimg = ifft2(ifftshift(circshift(FFimg, -[tX-1 tY-1])));

        Fimg = Fimg(3:end-2, 3:end-2);
        retAmplitude(:,:,iter) = single(abs(Fimg));

        % phase pipeline
        p = PhiShift(unwrap2(double(angle(Fimg))));
        p = phaseCompensation(p, 1);
        pp = p;
        tempThresh = p - imtophat(pp, diskSE);
        tempThresh = mean(tempThresh(:)) + 1.0;
        p2Mask = bwareaopen(pp > tempThresh(1), 100);
        p2Mask = ~imdilate(p2Mask, dilSE);
        p = phaseCompensation(p, 2, p2Mask);
        pNeg = p < 0;
        p = p - sum(p(pNeg)) / sum(pNeg(:));
        retPhase(:,:,iter) = single(p);
    end
    f_dx = f_dx_s;
    f_dy = f_dy_s;
    logfn(sprintf('  All %d frames processed.', nFrames));

    xSize = ii;
    [~, baseName, ~] = fileparts(sName);
    fileName = strcat('Field_', baseName);
    matOut = fullfile(outDir, strcat(fileName, '.mat'));
    pngOut  = fullfile(outDir, strcat(fileName, '.png'));
    logfn(sprintf('  Saving: %s', matOut));
    save(matOut, 'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
    logfn('  MAT file saved. Saving phase overview image...');
    saveFieldPNG(retPhase, nFrames, pngOut);
    logfn(sprintf('  Saved PNG: %s', pngOut));
end

%% Field inspection plot
logfn('Generating field inspection plot...');
fieldList = dir(fullfile(outDir, 'Field*.mat'));
logfn(sprintf('Found %d Field*.mat file(s) for inspection.', length(fieldList)));
inspectionPng = fullfile(outDir, 'Field_inspection.png');
saveInspectionPNG(outDir, fieldList, inspectionPng);
logfn(sprintf('Inspection plot saved: %s', inspectionPng));

logfn('=== field_Retrieval finished ===');
fclose(fid);
