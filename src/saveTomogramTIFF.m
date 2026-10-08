function saveTomogramTIFF(Reconimg, tifOut, flipLR)
% saveTomogramTIFF  Save Reconimg as a multi-page uint16 TIFF (scaled x10000).
%   flipLR (default true) mirrors each slice left-right before writing.
    if nargin < 3; flipLR = true; end
    try
        Reconimg_u16 = uint16(real(Reconimg) * 10000);
        if exist(tifOut, 'file'); delete(tifOut); end
        for K = 1:size(Reconimg_u16, 3)
            slice = Reconimg_u16(:,:,K);
            if flipLR; slice = fliplr(slice); end
            if K == 1
                imwrite(slice, tifOut, 'WriteMode', 'overwrite');
            else
                imwrite(slice, tifOut, 'WriteMode', 'append');
            end
        end
    catch ME
        warning('saveTomogramTIFF: %s', ME.message);
    end
end
