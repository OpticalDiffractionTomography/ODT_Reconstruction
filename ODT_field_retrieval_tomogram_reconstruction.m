%% Initial parameters
lambda = 0.532;     % [?m]  wavelength of the laser source
res = 4.8/57;
NA = 1.2;           % [1]   numerical apterature of the objective
numberImages = 150; % [1]   number of the acquired images

n_m = 1.3365;       % [1]   refractive index of the medium
n_s = 1.430;        % [1]   upper limit of the RI for the color bar (visualization)

folders = { ...     % folders to evaluate
    '.' ...
};

backgroundMeasurement = 'Background.h5';
measurements = 'Brillouin-3*';

for lll = 1:length(folders)
    folder = folders{lll};

    %% Background loading: c++ acquisition
    file = h5bmread([folder filesep 'RawData' filesep backgroundMeasurement]); %% Assign Background (no sample) fields
    % get the version attribute
    version = file.version;
    % get the date attribute
    date = file.date;

    ODTpayload = struct();
    ODTpayload.data = file.readPayloadData('ODT', 0, 'data', 48);
    tempData = ODTpayload.data;
    tomogMap = zeros([size(tempData), numberImages],'single');

    for kk = 1:numberImages
        ODTpayload.data = file.readPayloadData('ODT', 0, 'data',kk-1);
        tomogMap(:,:,kk) = ODTpayload.data;
    end
    %% Retrieve Background Field (no sample field)
    img = double(squeeze(tomogMap(:,:,49)));
    img = squeeze(img);
    ii = length(img);

    ZP = ii;
    r = round(ZP*res*NA/lambda) + 20;
    yr = r;
    Fimg = fftshift(fft2(img)) / (ii^2); %FFT
    [mi,mj] = find(Fimg==max(max(Fimg(:, round(ii*0.55):round(ii*1)))));
    mi = round(mi - ii/2-1);
    mj = round(mj - ii/2-1);

    c1mask = ~(mk_ellipse(r-20, yr-20, ZP, ZP));
    c3mask = circshift(c1mask, [mi mj]);

    Fbg = zeros(round(2*r), round(2*r), numberImages, 'single');
    tempImg = zeros(ZP,ZP);

    for bgNum = 1:numberImages
        img = double(squeeze(tomogMap(:,:,bgNum)));
        [xSize, ySize] = size(img);
        Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT    
        Fimg = Fimg.*c3mask;
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);

        Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
        sizeFimg = length(Fimg);
        Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
        Fbg(:,:,bgNum) = Fimg;
    end
    h5bmclose(file);
    %%
    sampleList = dir([folder filesep 'RawData' filesep measurements '.h5']);
    for sampleNum = 1:length(sampleList)
        try
            %% Field Loading: c++ acquisition
            [filePath, fileName, fileExt] = fileparts(sampleList(sampleNum).name);
            file = h5bmread([folder filesep 'RawData' filesep fileName fileExt]);
            repetitions = file.getRepetitions('ODT');
            for repetition = 1:length(repetitions)

                ODTpayload = struct();
                ODTpayload.data = file.readPayloadData('ODT', repetitions{repetition}, 'data', 1);
                tempData = ODTpayload.data;
                tomogMap = zeros([size(tempData), numberImages],'single');

                for kk = 1:numberImages
                    ODTpayload.data = file.readPayloadData('ODT', repetitions{repetition}, 'data', kk-1);
                    tomogMap(:,:,kk) = ODTpayload.data;
                end
                %%
                retPhase = zeros(round(2*r)-4, round(2*r)-4, numberImages, 'single');
                retAmplitude = zeros(round(2*r)-4, round(2*r)-4, numberImages, 'single');

                f_dx = NaN(numberImages, 1);
                f_dy = NaN(numberImages, 1);

                for iter = 1:numberImages
                    img = double(squeeze(tomogMap(:,:,iter)));
                    Fimg = fftshift(fft2(img))/(ii^2); %FFT            
                    [f_dx(iter), f_dy(iter)] = find(Fimg==max(max(Fimg(:, round(ii*0.55):round(ii*1) ))));
                    Fimg = Fimg.*c3mask;
                    Fimg = circshift(Fimg,-[round(mi) round(mj)]);
                    Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
                    sizeFimg = length(Fimg);            
                    Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
                    %% 
                    Fimg = Fimg ./ squeeze(Fbg(:,:,iter)); % Compensate by background field
                    Fimg = Fimg(3:end-2, 3:end-2);
                    retAmplitude(:,:,iter) = abs(Fimg);
                    p = PhiShift(unwrap2(double(angle(Fimg)))); % Compensate additional phase gradient/curvature
                    pmap = zeros(size(p));
                    p = (PhiShiftJH(p,~pmap));
                    p = p - median(median(p));
                    pmap = (abs(p)>1);
                    pmap = ~imdilate(pmap,strel('disk',5));
                    p = p - sum(sum(p.*pmap)) ./ sum(sum(pmap));
                    pmap = imgaussfilt(p .* pmap, numberImages);
                    ptemp = p-pmap;

                    level = multithresh(ptemp,1); % Make phase value outside sample
                    pmap = (ptemp) > level(1);
                    str = strel('disk',10);
                    pmap = imdilate(pmap,str);
                    ptemp = ptemp - sum(sum(ptemp.*~pmap)) ./ sum(sum(~pmap));
                    retPhase(:,:,iter) = ptemp;
                end

                if (exist('fig', 'var') && ishandle(fig))
                    figure(fig);
                else
                    fig = figure();
                end
                for kk = 1:6
                    subplot(3,3,kk);
                    imagesc(squeeze(retPhase(:,:,round(numberImages/6*kk))),[-3 3]);
                    axis image;
                    axis off;
                end
                subplot(3,3,[7 9]);
                imagesc(squeeze(retPhase(:,end/2,:)),[-3 3]);
                colormap('jet');

                fullFileName = strcat(folder, filesep, 'EvalData', filesep, 'Field_', fileName);
                if length(repetitions) > 1
                    fullFileName = [fullFileName '_rep' repetitions{repetition}];
                end
                save(strcat(fullFileName, '.mat'), 'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
                print('-dtiff', strcat(fullFileName, '.tif'));
            end
            h5bmclose(file);
        catch
        end
    end

    %%
    g = gpuDevice(1);
    reset(g);
    sampleList = dir([folder filesep 'EvalData' filesep 'Field_' measurements '.mat']);

    for sampleNum=1:length(sampleList)
        [filePath, fileName, fileExt] = fileparts(sampleList(sampleNum).name);
        load([folder filesep 'EvalData' filesep fileName fileExt]);

        [xx, yy, frame] = size(retPhase);
        crop_size = xx;

        f_dx2 = f_dx - f_dx(49); % subtract maxpoint
        f_dy2 = f_dy - f_dy(49);

        original_size=xSize;
        excludeFrame = [];
        temp = mean(squeeze(mean(abs(retPhase),1)));
        excludeFrame = [excludeFrame, find(abs(temp)>2 )]; 
        temp = temp - circshift(temp,1);
        excludeFrame = [excludeFrame, find(abs(temp)>0.05)];

        for kkk = 1:frame
            p2 = squeeze(retPhase(:,:,kkk));
            if isnan(max(max(p2)))
                excludeFrame = [excludeFrame, kkk];
            end
        end
        excludeFrame = unique(sort(excludeFrame));  % Exclude frames with abnormally high phase value or phase value differences due to dusts
        ZP = round(1.5*xx/2)*2; % Zero-padding in the E-field to increase
        crop_factor = crop_size/original_size;
        res2 = res/crop_factor;

        padd_factor = ZP/crop_size;
        kres = 1/(res*ZP)*crop_factor;
        f_dx2 = f_dx2*padd_factor;
        f_dy2 = f_dy2*padd_factor; 

        k0_x = kres*f_dx2; % for actual value, multiply resolution
        k0_y = kres*f_dy2;

        k0 = 1/lambda;
        k0_z = real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2)); % magnitude of absolute value is k0

        %%
        ZP2 = 512; % The size of 3-D Fourier space in the lateral direction
        ZP3 = 256; % The size of 3-D Fourier space in the axial direction
        res3 = res2*ZP/ZP2;
        res4 = res2*ZP/ZP3;

        ORytov = gpuArray(single(zeros(ZP2,ZP2,ZP3)));
        Count = (single(zeros(ZP2,ZP2,ZP3))); %Count is for averaging out multiple mappings at one pixel.
        frameList = 1:frame;
        frameList(excludeFrame) = [];
        for kk = frameList
            Hfilter = ones(ZP, ZP);
            FRytov = squeeze(log(retAmplitude(:,:,kk)) + 1i * retPhase(:,:,kk));
            FRytov = gpuArray(padarray(FRytov,[round((ZP-xx)/2) round((ZP-yy)/2)],'symmetric'));

            UsRytov = fftshift(fft2(FRytov)).*(res2)^2.*Hfilter;

            UsRytov = circshift(UsRytov,[round(f_dx2(kk)) round(f_dy2(kk))]);
            xr = ZP * res2 * NA / lambda;
            UsRytov = UsRytov .* ~mk_ellipse(xr, xr, ZP, ZP);
            [ky, kx] = meshgrid(kres*(-floor(ZP/2)+1:floor(ZP/2)),kres*(-floor(ZP/2)+1:floor(ZP/2)));
            kz = real(sqrt((n_m*k0)^2-kx.^2-ky.^2)); % Generating coordinates on the surface of Ewald Sphere 

            Kx = kx - k0_x(kk);
            Ky = ky - k0_y(kk);
            Kz = kz - k0_z(kk);
            Uprime = 1i .* 4 * pi * kz .* UsRytov; % Applying Fourier diffraction theorem

            xind = find((kz>0) .* ~mk_ellipse(xr,xr,ZP,ZP)...
                .*(Kx>(kres*(-floor(ZP2/2)+1)))...
                .*(Ky>(kres*(-floor(ZP2/2)+1)))...
                .*(Kz>(kres*(-floor(ZP3/2)+1)))...
                .*(Kx<(kres*(floor(ZP2/2))))...
                .*(Ky<(kres*(floor(ZP2/2))))...
                .*(Kz<(kres*(floor(ZP3/2))))); % Exclude information outside Ewald Sphere

            Uprime = Uprime(xind);
            Kx = Kx(xind);
            Ky = Ky(xind);
            Kz = Kz(xind);

            Kx = round(Kx/kres+ZP2/2);
            Ky = round(Ky/kres+ZP2/2);
            Kz = round(Kz/kres+ZP3/2);
            Kzp = (Kz-1)*ZP2^2+(Ky-1)*ZP2+Kx;
            temp = ORytov(Kzp);
            ORytov(Kzp) = temp+Uprime;
            Count(Kzp) = Count(Kzp)+1;
        end

        ORytov = gather(ORytov);
        ORytov(Count>0) = ORytov(Count>0)./Count(Count>0);
        ORytov = gpuArray(ORytov);
        clear Count UsRytov Uprime temp FRytov;
        %% Non-Negativity Constraint
        Reconimg = ifftn(ifftshift(ORytov))./(res3^2*res4);
        Reconimg = n_m*sqrt(1-Reconimg.*(lambda/(n_m*2*pi))^2);
        ORytov_index = (find(abs(ORytov)>0));
        ORytov = ORytov(ORytov_index);
        normFact = 1./(res3^2*res4);
        for mm = 1:100
            id = (real(Reconimg)<n_m);
            Reconimg(id) = n_m-1i*imag(Reconimg(id));
            clear id;
            Reconimg = -(2*pi*n_m/lambda)^2.*(Reconimg.^2/n_m^2-1);
            ORytov_new = fftshift(fftn(Reconimg))/normFact;
            ORytov_new(ORytov_index) = ORytov;
            Reconimg = (ifftn(ifftshift(ORytov_new)))*normFact;
            Reconimg = n_m*sqrt(1-Reconimg.*(lambda/(n_m*2*pi))^2);
        end
        Reconimg = real(gather(fftshift(Reconimg,3)));
        ZP4 = round(size(Reconimg,1)/3);
        ZP5 = round(size(Reconimg,3)/6);
        Reconimg = Reconimg(end/2-ZP4+1:end/2+ZP4,end/2-ZP4+1:end/2+ZP4,end/2-ZP5+1:end/2+ZP5);
        ZP4 = round(size(Reconimg,1));
        ZP5 = round(size(Reconimg,3));

        %% Plot
        if (exist('fig2', 'var') && ishandle(fig2))
            figure(fig2);
        else
            fig2 = figure();
        end
        subplot(221);
        imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-0,:,:)))'),[n_m-0.005 n_s]);
        axis image;
        colorbar;

        subplot(222);
        imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-0,:)))'),[n_m-0.005 n_s]);
        axis image;
        colorbar;

        subplot(223);
        imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+-1)))),[n_m-0.005 n_s]);
        axis image;
        colorbar;
        set(gca,'YDir','normal');

        subplot(224); % Maximum Projection of RI
        imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,max(real(Reconimg),[],3),[1.34 n_s])
        axis image;
        colorbar;
        set(gca,'YDir','normal');
        colormap('jet');

        %%
        fullFileName = strcat(folder, filesep, 'EvalData', filesep, 'Tomogram_', fileName);
        save(strcat(fullFileName, '.mat'),'Reconimg','res3','res4','lambda');
        print('-dtiff', strcat(fullFileName, '.tif'))
        close(fig2);
    end
end