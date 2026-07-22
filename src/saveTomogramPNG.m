function saveTomogramPNG(Reconimg, n_s, pngOut)
% saveTomogramPNG  Save a diagnostic PNG mosaic of tomogram slices.
    try
        clim  = [1.337-0.005, n_s];
        cmap  = jet(256);
        GAP   = 30;
        CBAR  = 24;
        ROWW  = 900;
        PH2   = 400;

        toRGB   = @(im) ind2rgb(im2uint8(mat2gray(im, clim)), cmap);
        resizeW = @(im, w) imresize(toRGB(im), [round(w*size(im,1)/size(im,2)), w]);

        xz  = real(squeeze(Reconimg(end/2,:,:)))';
        yz  = real(squeeze(Reconimg(:,end/2,:)))';
        xy  = real(squeeze(Reconimg(:,:,end/2+1)));
        mip = max(real(Reconimg), [], 3);

        panW1 = floor((ROWW - GAP) / 2);
        p_xz  = resizeW(xz, panW1);
        p_yz  = resizeW(yz, panW1);
        PH1   = max(size(p_xz,1), size(p_yz,1));
        wpad  = @(im, h) padarray(im, [h-size(im,1), 0], 1, 'post');
        p_xz  = wpad(p_xz, PH1);  p_yz = wpad(p_yz, PH1);
        hgap1 = ones(PH1, GAP, 3);
        row1  = [p_xz, hgap1, p_yz];

        p_xy  = imresize(toRGB(xy),  [PH2 PH2]);
        p_mip = imresize(toRGB(mip), [PH2 PH2]);
        hgap2 = ones(PH2, GAP, 3);
        row2_inner = [p_xy, hgap2, p_mip];
        pad2  = floor((ROWW - size(row2_inner,2)) / 2);
        pad2  = max(pad2, 0);
        row2  = [ones(PH2, pad2, 3), row2_inner, ones(PH2, ROWW-size(row2_inner,2)-pad2, 3)];

        vgap  = ones(GAP, ROWW, 3);
        grid  = [row1; vgap; row2];

        cbH   = size(grid, 1);
        cbar  = ind2rgb(im2uint8(repmat(linspace(1,0,cbH)', 1, CBAR)), cmap);
        cbgap = ones(cbH, GAP, 3);
        inner = [grid, cbgap, cbar];

        BORDER = 40;
        bH = size(inner,1); bW = size(inner,2);
        mosaic = ones(bH+2*BORDER, bW+2*BORDER, 3);
        mosaic(BORDER+1:BORDER+bH, BORDER+1:BORDER+bW, :) = inner;
        imwrite(mosaic, pngOut);
    catch ME
        warning('saveTomogramPNG: %s', ME.message);
    end
end
