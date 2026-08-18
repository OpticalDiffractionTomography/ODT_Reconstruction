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

%% Robustness limits for strongly scattering samples
% unwrap2 (Goldstein) can spin effectively forever on residue-heavy frames.
% Frames above residueLimit are routed to a serial watchdog pass and given at
% most unwrapTimeout seconds each; a sample exceeding sampleTimeLimit is
% abandoned so it cannot consume the whole SLURM walltime.
% Override any of these from main.sh / bash_template.sh via matlab -batch.
if ~exist('residueLimit', 'var');    residueLimit    = 20000; end  % residues/frame
if ~exist('unwrapTimeout', 'var');   unwrapTimeout   = 120;   end  % s per frame
if ~exist('sampleTimeLimit', 'var'); sampleTimeLimit = 3600;  end  % s per sample
logfn(sprintf('Limits: residueLimit=%d | unwrapTimeout=%g s | sampleTimeLimit=%g s', ...
    residueLimit, unwrapTimeout, sampleTimeLimit));
skippedSamples = {};

%% Start parallel pool (workers = SLURM CPUs minus 1 for main thread)
nWorkers = max(1, str2double(getenv('SLURM_CPUS_PER_TASK')) - 1);
if isnan(nWorkers) || nWorkers < 1, nWorkers = 4; end
if isempty(gcp('nocreate'))
    parpool('local', nWorkers);
end
logfn(sprintf('Parallel pool: %d workers', nWorkers));

%% File discovery
bglist = dir(fullfile(spath, 'bg*_Tomog.mat'));
logfn(sprintf('Found %d background file(s).', length(bglist)));
sampleList = dir(fullfile(spath, 'sample*_Tomog.mat'));
logfn(sprintf('Found %d sample file(s).', length(sampleList)));

% Every sample is paired with every background (matches original script);
% outputs are suffixed _bg_<N>.
for bgListNum = 1:length(bglist)
    bgName = bglist(bgListNum).name;
    logfn(sprintf('=== Background %d/%d: %s ===', bgListNum, length(bglist), bgName));
    load(fullfile(spath, bgName));
    img = double(squeeze(tomogMap(:,:,49)));
    img = squeeze(img);
    ii = length(img);
    ZP = ii;
    r = round(ZP*res*NA/lambda)+20;
    yr = r;
    logfn(sprintf('Image size: %d px | ZP=%d | r=%d | NA=%.3f | lambda=%.4f um | res=%.4f um', ii, ZP, r, NA, lambda, res));

    %% Carrier frequency detection from background frames
    nBg = size(tomogMap,3);
    logfn(sprintf('Detecting carrier frequency from %d background frames...', nBg));
    f_dx = zeros(nBg,1);
    f_dy = zeros(nBg,1);
    parfor bgNum = 1:nBg
        bgImg = double(squeeze(tomogMap(:,:,bgNum)));
        Ftmp = fftshift(fft2(bgImg))/(ii^2);
        [fx,fy] = find(Ftmp==max(max(Ftmp(:, round(ii*0.01):round(ii*0.49)))));
        f_dx(bgNum) = fx(1);
        f_dy(bgNum) = fy(1);
    end
    % Carrier offset taken from frame 49 (matches original script)
    mi = f_dx(49);
    mj = f_dy(49);
    mi = round(mi-ii/2-1); mj = round(mj-ii/2-1);
    logfn(sprintf('Carrier frequency offset (frame 49): mi=%d px, mj=%d px', mi, mj));

    %% Build demodulation mask and background field stack
    logfn('Building demodulation mask (mk_ellipse)...');
    c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
    c3mask = circshift(c1mask,[mi mj]);
    logfn(sprintf('Demodulating %d background frames...', nBg));
    Fbg = zeros(round(2*r),round(2*r),nBg,'single');
    mi_r = round(mi); mj_r = round(mj);
    rowRange = ii/2-r+1:ii/2+r;
    parfor bgNum = 1:nBg
        bgImg = double(squeeze(tomogMap(:,:,bgNum)));
        Ftmp = fftshift(fft2(bgImg))/(ii^2);
        Ftmp = Ftmp.*c3mask;
        Ftmp = circshift(Ftmp,-[mi_r mj_r]);
        Ftmp = Ftmp(rowRange, rowRange);
        sz   = length(Ftmp);
        Fbg(:,:,bgNum) = single(ifft2(ifftshift(Ftmp))*(sz^2));
    end
    logfn('Background field stack ready.');

    %% Per-sample field retrieval
    for sampleNum = 1:length(sampleList)
        sName = sampleList(sampleNum).name;
        logfn(sprintf('--- Processing sample %d/%d (bg %d): %s ---', sampleNum, length(sampleList), bgListNum, sName));

        tSample = tic;
        logfn(sprintf('  Loading %s...', sName));
        load(fullfile(spath, sName));
        nFrames = size(tomogMap,3);
        logfn(sprintf('  Loaded: %d frames, size %dx%d', nFrames, size(tomogMap,1), size(tomogMap,2)));

        pSize = round(2*r)-4;
        retPhase     = zeros(pSize, pSize, nFrames, 'single');
        retAmplitude = zeros(pSize, pSize, nFrames, 'single');
        f_dx_s = zeros(nFrames,1);
        f_dy_s = zeros(nFrames,1);

        ii_sq  = ii^2;
        colLo  = round(ii*0.01);
        colHi  = round(ii*0.49);
        diskSE = strel('disk', 150);  % created once; safe to broadcast (value object)
        dilSE  = strel('disk', 5);

        % Pass 1 (parallel): demodulation + phase pipeline for every frame.
        % A frame whose unwrap2 call hangs cannot be interrupted from inside
        % parfor, so instead of calling unwrap2 here for frames that look
        % dangerous, the wrapped phase is stashed and retried serially in
        % pass 2 under a wall-clock watchdog.
        % Deferred frames wrapped phase is kept in a sliced cell (only the
        % few risky frames occupy memory, not a second full-size stack).
        demodPhase   = cell(nFrames,1);
        deferFrame   = false(nFrames,1);
        residueCount = zeros(nFrames,1);
        logfn(sprintf('  Launching parfor over %d frames (%d workers)...', nFrames, nWorkers));
        parfor iter = 1:nFrames
            frImg = double(squeeze(tomogMap(:,:,iter)));
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
            Fimg = Fimg ./ squeeze(Fbg(:,:,iter));

            % residual tilt correction via peak shift
            FFimg = fftshift(fft2(Fimg));
            [tX,tY] = find(abs(FFimg) == max(max(abs(FFimg))));
            tX = tX(1) - size(FFimg,1)/2;
            tY = tY(1) - size(FFimg,2)/2;
            Fimg = ifft2(ifftshift(circshift(FFimg, -[tX-1 tY-1])));

            Fimg = Fimg(3:end-2, 3:end-2);
            retAmplitude(:,:,iter) = single(abs(Fimg));

            wrappedPhase = double(angle(Fimg));

            % Residue density predicts unwrap2 blow-ups: Goldstein pairs
            % residues into branch cuts, and its cost explodes when there are
            % many. Strongly scattering samples produce exactly that.
            nRes = phaseResidueCount(wrappedPhase);
            residueCount(iter) = nRes;

            if any(~isfinite(wrappedPhase(:))) || nRes > residueLimit
                % Defer: unwrap2 either cannot run or is likely to hang.
                deferFrame(iter) = true;
                demodPhase{iter} = single(wrappedPhase);
                retPhase(:,:,iter) = NaN(pSize, pSize, 'single');
            else
                try
                    retPhase(:,:,iter) = single(phasePipeline(wrappedPhase, diskSE, dilSE));
                catch ME
                    fprintf('phase pipeline failed at iter = %d: %s\n', iter, ME.message);
                    deferFrame(iter) = true;
                    demodPhase{iter} = single(wrappedPhase);
                    retPhase(:,:,iter) = NaN(pSize, pSize, 'single');
                end
            end
        end
        f_dx = f_dx_s;
        f_dy = f_dy_s;
        logfn(sprintf('  Pass 1 done. Residues per frame: median %d, max %d (limit %d).', ...
            round(median(residueCount)), max(residueCount), residueLimit));

        % Pass 2 (serial, watchdogged): retry the deferred frames one at a
        % time. parfeval + wait() gives a hard time limit that a try/catch
        % cannot, so a hanging unwrap2 costs `unwrapTimeout` seconds instead
        % of stalling the whole job.
        deferList = find(deferFrame(:)');
        if ~isempty(deferList)
            logfn(sprintf('  Pass 2: retrying %d deferred frame(s) with %g s timeout each: %s', ...
                numel(deferList), unwrapTimeout, num2str(deferList)));
            for iter = deferList
                if toc(tSample) > sampleTimeLimit
                    logfn(sprintf('    sample time limit (%g s) reached — abandoning retries for remaining frames', ...
                        sampleTimeLimit));
                    skippedSamples{end+1} = sName;  %#ok<SAGROW>
                    break
                end
                wrappedPhase = double(demodPhase{iter});
                if isempty(wrappedPhase) || any(~isfinite(wrappedPhase(:)))
                    logfn(sprintf('    frame %d: non-finite wrapped phase — skipped', iter));
                    continue
                end
                tRetry = tic;
                try
                    p = unwrapWithTimeout(wrappedPhase, unwrapTimeout);
                    retPhase(:,:,iter) = single(phasePipeline(p, diskSE, dilSE, true));
                    logfn(sprintf('    frame %d: recovered in %.1f s (%d residues)', ...
                        iter, toc(tRetry), residueCount(iter)));
                catch ME
                    logfn(sprintf('    frame %d: %s after %.1f s (%d residues) — NaN', ...
                        iter, ME.identifier, toc(tRetry), residueCount(iter)));
                end
            end
        end

        badFrames = find(squeeze(any(any(isnan(retPhase), 1), 2)));
        if ~isempty(badFrames)
            logfn(sprintf('  WARNING: %d/%d frame(s) unusable (NaN, excluded downstream): %s', ...
                numel(badFrames), nFrames, num2str(badFrames(:)')));
        end
        if numel(badFrames) > 0.5*nFrames
            logfn(sprintf('  WARNING: >50%% of frames unusable — %s is likely too strongly scattering to reconstruct reliably.', sName));
        end
        logfn(sprintf('  All %d frames processed.', nFrames));

        xSize = ii;
        [~, baseName, ~] = fileparts(sName);
        fileName = strcat('Field_', baseName, '_bg_', num2str(bgListNum));
        matOut = fullfile(outDir, strcat(fileName, '.mat'));
        pngOut  = fullfile(outDir, strcat(fileName, '.png'));
        logfn(sprintf('  Saving: %s', matOut));
        save(matOut, 'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
        logfn('  MAT file saved. Saving phase overview image...');
        saveFieldPNG(retPhase, nFrames, pngOut);
        logfn(sprintf('  Saved PNG: %s', pngOut));
    end
end

%% Field inspection plot
logfn('Generating field inspection plot...');
fieldList = dir(fullfile(outDir, 'Field*.mat'));
logfn(sprintf('Found %d Field*.mat file(s) for inspection.', length(fieldList)));
inspectionPng = fullfile(outDir, 'Field_inspection.png');
saveInspectionPNG(outDir, fieldList, inspectionPng);
logfn(sprintf('Inspection plot saved: %s', inspectionPng));

if ~isempty(skippedSamples)
    logfn(sprintf('PROBLEM SAMPLES (hit the %g s limit; some frames NaN): %s', ...
        sampleTimeLimit, strjoin(unique(skippedSamples), ', ')));
end

logfn('=== field_Retrieval finished ===');
fclose(fid);
