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
    if ~exist('n_m', 'var'); n_m = 1.337; end
    n_s = n_m + 0.04;
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
    logfn('  MAT file saved.');

    tifOut = fullfile(outDir, strcat(fileName, '.tif'));
    logfn(sprintf('  Saving TIFF: %s', tifOut));
    saveTomogramTIFF(Reconimg, tifOut);
    logfn(sprintf('  TIFF saved (%d slices).', size(Reconimg,3)));

    logfn(sprintf('  Saving PNG: %s', pngOut));
    saveTomogramPNG(Reconimg, n_s, pngOut);
    logfn('  PNG saved.');
end

logfn('=== tomogram_Reconstruction finished ===');
fclose(fid);
