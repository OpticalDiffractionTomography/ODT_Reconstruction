function saveFieldPNG(retPhase, nFrames, pngOut)
% saveFieldPNG  Save a diagnostic PNG: 6 phase panels + kymograph strip.
    try
        cmap    = jet(256);
        GAP     = 20;
        CBAR    = 20;
        P       = 300;
        clim_ph = [-3 3];
        toRGB   = @(im) ind2rgb(im2uint8(mat2gray(im, clim_ph)), cmap);

        panels = cell(1,6);
        for kk = 1:6
            panels{kk} = imresize(toRGB(squeeze(retPhase(:,:,round(nFrames/6*kk)))), [P P]);
        end

        rowW   = 3*P + 2*GAP;
        Wtotal = rowW + GAP + CBAR;
        hgap   = ones(P, GAP, 3);
        row1   = [panels{1}, hgap, panels{2}, hgap, panels{3}];
        row2   = [panels{4}, hgap, panels{5}, hgap, panels{6}];
        row1   = padarray(row1, [0, Wtotal-rowW], 1, 'post');
        row2   = padarray(row2, [0, Wtotal-rowW], 1, 'post');

        kymoH     = 200;
        strip_rs  = imresize(toRGB(squeeze(retPhase(:,end/2,:))), [kymoH, rowW]);
        cbar_st   = ind2rgb(im2uint8(repmat(linspace(1,0,kymoH)', 1, CBAR)), cmap);
        strip_row = [strip_rs, ones(kymoH, GAP, 3), cbar_st];

        vgap   = ones(GAP, Wtotal, 3);
        grid   = [row1; vgap; row2; vgap; strip_row];

        BORDER = 40;
        bH = size(grid,1); bW = size(grid,2);
        mosaic = ones(bH+2*BORDER, bW+2*BORDER, 3);
        mosaic(BORDER+1:BORDER+bH, BORDER+1:BORDER+bW, :) = grid;
        imwrite(mosaic, pngOut);
    catch ME
        warning('saveFieldPNG: %s', ME.message);
    end
end
