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
for bgNum = 1:size(tomogMap,3)
    img = double(squeeze(tomogMap(:,:,bgNum)));
    [xSize ySize] = size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize);
    [f_x,f_y] = find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49)))));
    f_dx = [f_dx; f_x];
    f_dy = [f_dy; f_y];
end
mi = mean(f_dx);
mj = mean(f_dy);
mi = round(mi-ii/2-1); mj = round(mj-ii/2-1);
logfn(sprintf('Carrier frequency offset: mi=%d px, mj=%d px', mi, mj));

%% Build demodulation mask and background field stack
logfn('Building demodulation mask (mk_ellipse)...');
c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
c3mask = circshift(c1mask,[mi mj]);
logfn(sprintf('Demodulating %d background frames...', size(tomogMap,3)));
Fbg = zeros(round(2*r),round(2*r),size(tomogMap,3),'single');
for bgNum = 1:size(tomogMap,3)
    img = double(squeeze(tomogMap(:,:,bgNum)));
    [xSize ySize] = size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize);
    Fimg = Fimg.*c3mask;
    Fimg = circshift(Fimg,-[round(mi) round(mj)]);
    Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
    sizeFimg = length(Fimg);
    Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
    Fbg(:,:,bgNum) = Fimg;
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

    retPhase     = zeros(round(2*r)-4, round(2*r)-4, nFrames, 'single');
    retAmplitude = zeros(round(2*r)-4, round(2*r)-4, nFrames, 'single');
    f_dx = [];
    f_dy = [];

    for iter = 1:nFrames
        if mod(iter,10)==1 || iter==nFrames
            logfn(sprintf('  Frame %d/%d: demodulation + phase retrieval', iter, nFrames));
        end
        img = (double(squeeze(tomogMap(:,:,iter))));
        Fimg = fftshift(fft2(img))/(ii^2);
        [f_x,f_y] = find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49)))));
        f_dx = [f_dx; f_x];
        f_dy = [f_dy; f_y];
        Fimg = Fimg.*c3mask;
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);
        Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
        sizeFimg = length(Fimg);
        Fimg = (ifft2(ifftshift(Fimg))*(sizeFimg^2));
        Fimg = Fimg./squeeze(Fbg(:,:,iter));

        FFimg = (fftshift(fft2(Fimg)));
        [tempX,tempY] = find(abs(FFimg)==max(max(abs(FFimg))));
        tempX = tempX-size(FFimg,1)/2;
        tempY = tempY-size(FFimg,2)/2;
        Fimg = ifft2(ifftshift(circshift(FFimg,-[tempX-1 tempY-1])));

        Fimg = Fimg(3:end-2,3:end-2);
        retAmplitude(:,:,iter) = abs(Fimg);
        p = (PhiShift(unwrap2(double(angle(Fimg)))));
        [p,coeffVal,~] = phaseCompensation(p,1);
        pp = (p);
        tempThresh = p-imtophat(pp,strel('disk',150));
        tempThresh = mean(tempThresh(:))+1.0;
        p2Mask = (pp>tempThresh(1));
        p2Mask = bwareaopen(p2Mask,100);
        p2Mask = ~imdilate(p2Mask,strel('disk',5));
        [p,coefficients] = phaseCompensation(p,2,p2Mask);
        pNegative = (p<0);
        p = p-sum(sum(p(pNegative)))./sum(sum(pNegative));
        retPhase(:,:,iter) = p;
    end
    logfn(sprintf('  All %d frames processed.', nFrames));

    [~, baseName, ~] = fileparts(sName);
    fileName = strcat('Field_', baseName);
    matOut = fullfile(outDir, strcat(fileName, '.mat'));
    pngOut  = fullfile(outDir, strcat(fileName, '.png'));
    logfn(sprintf('  Saving: %s', matOut));
    save(matOut, 'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
    logfn('  MAT file saved. Saving phase overview image...');
    try
        cmap    = jet(256);
        GAP     = 20;       % white gap between panels (px)
        clim_ph = [-3 3];
        toRGB   = @(im) ind2rgb(im2uint8(mat2gray(im, clim_ph)), cmap);

        % 6 square phase panels — resize all to same size
        P = 300;   % each panel size in pixels
        panels = cell(1,6);
        for kk = 1:6
            panels{kk} = imresize(toRGB(squeeze(retPhase(:,:,round(nFrames/6*kk)))), [P P]);
        end

        CBAR  = 20;
        % total width = 3 panels + 2 gaps + gap + colorbar
        rowW  = 3*P + 2*GAP;
        Wtotal = rowW + GAP + CBAR;

        hgap  = ones(P, GAP, 3);
        row1  = [panels{1}, hgap, panels{2}, hgap, panels{3}];
        row2  = [panels{4}, hgap, panels{5}, hgap, panels{6}];
        % pad rows to Wtotal
        row1  = padarray(row1, [0, Wtotal-rowW], 1, 'post');
        row2  = padarray(row2, [0, Wtotal-rowW], 1, 'post');

        % kymograph strip: fixed height, stretched to full rowW width
        strip     = toRGB(squeeze(retPhase(:,end/2,:)));
        kymoH     = 200;
        strip_rs  = imresize(strip, [kymoH, rowW]);
        cbgap_st  = ones(kymoH, GAP, 3);
        cbar_st   = ind2rgb(im2uint8(repmat(linspace(1,0,kymoH)', 1, CBAR)), cmap);
        strip_row = [strip_rs, cbgap_st, cbar_st];  % width == Wtotal

        vgap = ones(GAP, Wtotal, 3);
        grid = [row1; vgap; row2; vgap; strip_row];

        % outer white border
        BORDER = 40;
        bH = size(grid,1); bW = size(grid,2);
        mosaic = ones(bH+2*BORDER, bW+2*BORDER, 3);
        mosaic(BORDER+1:BORDER+bH, BORDER+1:BORDER+bW, :) = grid;
        imwrite(mosaic, pngOut);
        logfn(sprintf('  Saved PNG: %s', pngOut));
    catch ME
        logfn(sprintf('  WARNING: PNG save failed (%s) — continuing.', ME.message));
    end
end

%% Field inspection plot (headless-safe: renders to pixel buffer, no Qt)
logfn('Generating field inspection plot...');
fieldList = dir(fullfile(outDir, 'Field*.mat'));
logfn(sprintf('Found %d Field*.mat file(s) for inspection.', length(fieldList)));
inspectionPng = fullfile(outDir, 'Field_inspection.png');
try
    W = 800; H = 400;
    canvas = ones(H, W, 3, 'uint8') * 255;  % white background
    allTemp = {}; allDiff = {};
    for sampleNum = 1:length(fieldList)
        load(fullfile(outDir, fieldList(sampleNum).name), 'retPhase');
        t = mean(squeeze(mean(abs(retPhase), 1)));
        allTemp{end+1} = t;
        allDiff{end+1} = t - circshift(t, 1);
    end
    ymin = -3; ymax = 3;
    for sampleNum = 1:length(allTemp)
        t = allTemp{sampleNum};
        d = allDiff{sampleNum};
        N = length(t);
        for fr = 1:N
            px = round((fr-1)/(N-1) * (W-1)) + 1;
            py_t = round((1 - (t(fr)-ymin)/(ymax-ymin)) * (H-1)) + 1;
            py_d = round((1 - (d(fr)-ymin)/(ymax-ymin)) * (H-1)) + 1;
            py_t = max(1, min(H, py_t));
            py_d = max(1, min(H, py_d));
            canvas(py_t, px, :) = reshape([255 0 0], 1, 1, 3);  % red: mean phase
            canvas(py_d, px, :) = reshape([0 180 0], 1, 1, 3);  % green: diff
        end
    end
    imwrite(canvas, inspectionPng);
    logfn(sprintf('Inspection plot saved: %s', inspectionPng));
catch ME
    logfn(sprintf('WARNING: inspection plot failed (%s) — continuing.', ME.message));
end

logfn('=== field_Retrieval finished ===');
fclose(fid);
