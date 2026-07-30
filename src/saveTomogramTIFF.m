function saveTomogramTIFF(Reconimg, tifOut)
% saveTomogramTIFF  Save Reconimg as a multi-page uint16 TIFF (scaled x10000).
    try
        Reconimg_u16 = uint16(real(Reconimg) * 10000);
        if exist(tifOut, 'file'); delete(tifOut); end
        for K = 1:size(Reconimg_u16, 3)
            if K == 1
                imwrite(fliplr(Reconimg_u16(:,:,K)), tifOut, 'WriteMode', 'overwrite');
            else
                imwrite(fliplr(Reconimg_u16(:,:,K)), tifOut, 'WriteMode', 'append');
            end
        end
    catch ME
        warning('saveTomogramTIFF: %s', ME.message);
    end
end
