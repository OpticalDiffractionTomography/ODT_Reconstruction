function saveInspectionPNG(outDir, fieldList, inspectionPng)
% saveInspectionPNG  Save mean-phase and frame-diff curves as a pixel-rendered PNG.
    try
        W = 800; H = 400;
        canvas = ones(H, W, 3, 'uint8') * 255;
        allTemp = {}; allDiff = {};
        for s = 1:length(fieldList)
            load(fullfile(outDir, fieldList(s).name), 'retPhase');
            t = mean(squeeze(mean(abs(retPhase), 1)));
            allTemp{end+1} = t;
            allDiff{end+1} = t - circshift(t, 1);
        end
        ymin = -3; ymax = 3;
        for s = 1:length(allTemp)
            t = allTemp{s};  d = allDiff{s};  N = length(t);
            for fr = 1:N
                px   = round((fr-1)/(N-1) * (W-1)) + 1;
                py_t = max(1, min(H, round((1-(t(fr)-ymin)/(ymax-ymin))*(H-1))+1));
                py_d = max(1, min(H, round((1-(d(fr)-ymin)/(ymax-ymin))*(H-1))+1));
                canvas(py_t, px, :) = reshape([255 0 0], 1, 1, 3);
                canvas(py_d, px, :) = reshape([0 180 0], 1, 1, 3);
            end
        end
        imwrite(canvas, inspectionPng);
    catch ME
        warning('saveInspectionPNG: %s', ME.message);
    end
end
