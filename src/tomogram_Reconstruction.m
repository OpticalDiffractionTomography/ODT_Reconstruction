%% Tomogram Reconstruction
% spath is provided by main.sh at invocation
set(0, 'DefaultFigureVisible', 'off');  % headless: no display on compute nodes

outDir = fullfile(spath, 'field_retrieval');

logFile = fullfile(outDir, 'tomogram_reconstruction.log');
fid = fopen(logFile, 'w');
logfn = @(msg) fprintf(fid, '[%s] %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'), msg);
logfn('=== tomogram_Reconstruction started ===');
logfn(sprintf('spath: %s', spath));

sampleList = dir(fullfile(outDir, 'Field*.mat'));
logfn(sprintf('Found %d Field*.mat file(s).', length(sampleList)));

for sampleNum = 1:length(sampleList)
    sName = sampleList(sampleNum).name;
    logfn(sprintf('--- Processing %d/%d: %s ---', sampleNum, length(sampleList), sName));

    logfn(sprintf('  Loading %s...', sName));
    load(fullfile(outDir, sName));
    [xx, yy, frame] = size(retPhase);
    logfn(sprintf('  Loaded: retPhase size %dx%dx%d', xx, yy, frame));

    crop_size = xx;
    f_dx2 = f_dx-mean(f_dx(:));
    f_dy2 = f_dy-mean(f_dy(:));
    original_size = xSize;

    %% Outlier frame detection
    logfn('  Detecting outlier frames...');
    excludeFrame = [];
    temp = mean(squeeze(mean(abs(retPhase),1)));
    excludeFrame = [excludeFrame, find(abs(temp)>1.5)];
    temp = temp-circshift(temp,1);
    excludeFrame = [excludeFrame, find(abs(temp)>0.1)];
    for kkk = 1:frame
        p2 = squeeze(retPhase(:,:,kkk));
        if sum(isnan(p2(:)))
            excludeFrame = [excludeFrame, kkk];
        end
    end
    excludeFrame = unique(sort(excludeFrame));
    logfn(sprintf('  Excluded %d frame(s): [%s]', length(excludeFrame), num2str(excludeFrame)));

    %% Build TomoParam
    logfn('  Building TomoParam struct...');
    n_m = 1.337;
    n_s = 1.337+0.04;
    ZP = round(1.2*xx/2)*2;
    crop_factor = crop_size/original_size;
    res2 = res/crop_factor;
    padd_factor = ZP/crop_size;
    kres = 1/(res*ZP)*crop_factor;
    f_dx2 = f_dx2*padd_factor; f_dy2 = f_dy2*padd_factor;
    k0_x = kres*f_dx2;
    k0_y = kres*f_dy2;
    k0 = 1/lambda;
    k0_z = real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2));
    ZP2 = round(512*1.0/2)*2;
    ZP3 = round(256*1.0/2)*2;
    res3 = res2*ZP/ZP2;
    res4 = res2*ZP/ZP3;
    frameList = (1:frame);
    frameList(excludeFrame) = [];
    IterNum = 100;
    logfn(sprintf('  ZP=%d | ZP2=%d | ZP3=%d | frameList length=%d', ZP, ZP2, ZP3, length(frameList)));
    logfn(sprintf('  res2=%.4f um | res3=%.4f um | res4=%.4f um', res2, res3, res4));
    logfn(sprintf('  n_m=%.4f | lambda=%.4f um | NA=%.3f', n_m, lambda, NA));

    TomoParam = struct('n_m',n_m,'ZP',ZP,'ZP2',ZP2,'ZP3',ZP3,'xx',xx,'yy',yy,...
        'f_dx2',f_dx2,'f_dy2',f_dy2,'NA',NA,'lambda',lambda,...
        'k0',k0,'k0_x',k0_x,'k0_y',k0_y,'k0_z',k0_z,'kres',kres,'frameList',frameList,...
        'res2',res2,'res3',res3,'res4',res4);

    %% Ewald sphere mapping
    logfn('  Running ODTReconstruction (Ewald sphere mapping)...');
    [Reconimg, ORytov] = ODTReconstruction(retAmplitude,retPhase,TomoParam);
    logfn('  ODTReconstruction done.');

    %% Iterative refinement
    logfn(sprintf('  Running ODTIteration (%d iterations on GPU)...', IterNum));
    Reconimg = ODTIteration(gpuArray((Reconimg)),ORytov,TomoParam,IterNum);
    logfn('  ODTIteration done.');

    %% Crop and save
    logfn('  Cropping to unpadded region...');
    ZP4 = round(size(Reconimg,1)/(2*padd_factor));
    ZP5 = round(size(Reconimg,3)/(4*padd_factor));
    Reconimg = real(Reconimg(end/2-ZP4+1:end/2+ZP4,end/2-ZP4+1:end/2+ZP4,end/2-ZP5+1:end/2+ZP5));
    ZP4 = round(size(Reconimg,1));
    ZP5 = round(size(Reconimg,3));
    logfn(sprintf('  Final Reconimg size: %dx%dx%d', ZP4, ZP4, ZP5));

    [~, baseName, ~] = fileparts(sName);
    fileName = strcat('Tomogram_', baseName);
    matOut = fullfile(outDir, strcat(fileName, '.mat'));
    pngOut  = fullfile(outDir, strcat(fileName, '.png'));
    logfn(sprintf('  Saving: %s', matOut));
    save(matOut, 'Reconimg','res3','res4','lambda','excludeFrame');
    logfn('  MAT file saved. Saving reconstruction image...');
    try
        clim  = [1.337-0.005, n_s];
        cmap  = jet(256);
        GAP   = 30;     % gap between panels and border (px)
        CBAR  = 24;     % colorbar width (px)
        ROWW  = 900;    % target total width per row (both panels + gap)
        PH2   = 400;    % row-2 panel height (XY / MIP — larger)

        toRGB   = @(im) ind2rgb(im2uint8(mat2gray(im, clim)), cmap);
        % resize panel to fixed height, preserving aspect ratio
        resizeH = @(im, h) imresize(toRGB(im), [h, round(h*size(im,2)/size(im,1))]);
        % resize panel to fixed width
        resizeW = @(im, w) imresize(toRGB(im), [round(w*size(im,1)/size(im,2)), w]);

        xz  = real(squeeze(Reconimg(end/2,:,:)))';
        yz  = real(squeeze(Reconimg(:,end/2,:)))';
        xy  = real(squeeze(Reconimg(:,:,end/2+1)));
        mip = max(real(Reconimg), [], 3);

        % Row 1: XZ and YZ — each gets half of ROWW (minus gap)
        panW1 = floor((ROWW - GAP) / 2);
        p_xz  = resizeW(xz, panW1);
        p_yz  = resizeW(yz, panW1);
        PH1   = max(size(p_xz,1), size(p_yz,1));
        wpad  = @(im, h) padarray(im, [h-size(im,1), 0], 1, 'post');
        p_xz  = wpad(p_xz, PH1);  p_yz = wpad(p_yz, PH1);
        hgap1 = ones(PH1, GAP, 3);
        row1  = [p_xz, hgap1, p_yz];

        % Row 2: XY and MIP — square panels at PH2 height, centered in ROWW
        p_xy  = imresize(toRGB(xy),  [PH2 PH2]);
        p_mip = imresize(toRGB(mip), [PH2 PH2]);
        hgap2 = ones(PH2, GAP, 3);
        row2_inner = [p_xy, hgap2, p_mip];
        % center row2 within ROWW
        pad2  = floor((ROWW - size(row2_inner,2)) / 2);
        pad2  = max(pad2, 0);
        row2  = [ones(PH2, pad2, 3), row2_inner, ones(PH2, ROWW-size(row2_inner,2)-pad2, 3)];

        vgap  = ones(GAP, ROWW, 3);
        grid  = [row1; vgap; row2];

        % single colorbar on the right spanning full grid height
        cbH   = size(grid, 1);
        cbar  = ind2rgb(im2uint8(repmat(linspace(1,0,cbH)', 1, CBAR)), cmap);
        cbgap = ones(cbH, GAP, 3);
        inner = [grid, cbgap, cbar];

        BORDER = 40;
        bH = size(inner,1); bW = size(inner,2);
        mosaic = ones(bH+2*BORDER, bW+2*BORDER, 3);
        mosaic(BORDER+1:BORDER+bH, BORDER+1:BORDER+bW, :) = inner;
        imwrite(mosaic, pngOut);
        logfn(sprintf('  Saved PNG: %s', pngOut));
    catch ME
        logfn(sprintf('  WARNING: PNG save failed (%s) — continuing.', ME.message));
    end
end

logfn('=== tomogram_Reconstruction finished ===');
fclose(fid);
