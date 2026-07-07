clear all
close all
addpath(genpath('/beegfs/home/ralajan/matlab'))
%% Angle load
spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK';
cd(spath);
batchList=dir('batch*');
for batchNum=1:length(batchList)
    cd(spath)
    cd(batchList(batchNum).name)
    f_dx =[];
    f_dy= [];
    bglist=dir('bg*_Tomog.mat');
    cd(spath)
    cd(batchList(batchNum).name)
    load(bglist(1).name);
    img=double(squeeze(tomogMap(:,:,49)));
    %                 img=double(imread(bglist(1).name,49));
    img=squeeze(img);
    ii=length(img);
    ZP=ii;
    r=round(ZP*res*NA/lambda)+20;
    yr = r;
    for bgNum=1:size(tomogMap,3)
        img=double(squeeze(tomogMap(:,:,bgNum)));
        %                 subplot(121),  imagesc(((img))),axis image,colormap('jet')
        [xSize ySize]=size(img);
        Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT
        %             subplot(122),imagesc(log10(abs(Fimg))),axis image
        [f_x,f_y]=find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49) ))));
        f_dx=[f_dx;f_x];
        f_dy=[f_dy;f_y];
        %             pause()
    end
    mi=f_dx(49);
    mj=f_dy(49);
    %         mi=mean(f_dx);
    %         mj=mean(f_dy);
    mi=round(mi-ii/2-1); mj=round(mj-ii/2-1);
    c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
    c3mask = circshift(c1mask,[mi mj]);
    Fbg=zeros(round(2*r),round(2*r),size(tomogMap,3),'single');
    tempImg=zeros(ZP,ZP);
    for bgNum=1:size(tomogMap,3)
        img=double(squeeze(tomogMap(:,:,bgNum)));
        [xSize ySize]=size(img);
        Fimg = fftshift(fft2(img))/(xSize*ySize); %FFT
        Fimg = Fimg.*c3mask;
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);
        Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
        sizeFimg=length(Fimg);
        Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
        Fbg(:,:,bgNum)=Fimg;
    end
    %%
    sampleList=dir('sample*_Tomog.mat');
    for sampleNum=1:length(sampleList)
        cd(spath)
        cd(batchList(batchNum).name)
        load(sampleList(sampleNum).name);
        retPhase=zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
        retAmplitude=zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
        f_dx =[];
        f_dy= [];
        for iter=1:size(tomogMap,3)
            img=(double(squeeze(tomogMap(:,:,iter))));
            %                     subplot(231),  imagesc(((img))),axis image
            Fimg = fftshift(fft2(img))/(ii^2); %FFT
            [f_x,f_y]=find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49) ))));
            f_dx=[f_dx;f_x];
            f_dy=[f_dy;f_y];
            %                     subplot(232),  imagesc(log10(abs(Fimg))),axis image
            %                     subplot(2,3,[4 6]),plot(iter,std(img(:)),'or'),hold on
            Fimg = Fimg.*c3mask;
            Fimg = circshift(Fimg,-[round(mi) round(mj)]);
            Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
            sizeFimg=length(Fimg);
            Fimg = (ifft2(ifftshift(Fimg))*(sizeFimg^2));
            Fimg=Fimg./squeeze(Fbg(:,:,iter));
            %%
            %                     figure,imagesc(log10(),axis image
            FFimg=(fftshift(fft2(Fimg)));
            [tempX,tempY]=find(abs(FFimg)==max(max(abs(FFimg))));
            tempX=tempX-size(FFimg,1)/2;
            tempY=tempY-size(FFimg,2)/2;
            Fimg=ifft2(ifftshift(circshift(FFimg,-[tempX-1 tempY-1])));
            %%
            Fimg=Fimg(3:end-2,3:end-2);
            retAmplitude(:,:,iter)=abs(Fimg);
            p=(PhiShift(unwrap2(double(angle(Fimg)))));
            [p,coeffVal,~]=phaseCompensation(p,1);
            %                     subplot(233),imagesc(p),axis image
            pp=(p);
            tempThresh=p-imtophat(pp,strel('disk',150));
            tempThresh=mean(tempThresh(:))+1.0;
            p2Mask=(pp>tempThresh(1));
            p2Mask=bwareaopen(p2Mask,100);
            p2Mask=~imdilate(p2Mask,strel('disk',5));
            [p,coefficients]=phaseCompensation(p,2,p2Mask);
            pNegative=(p<0);
            %             pNegative=imerode(pNegative,strel('disk',3));
            p=p-sum(sum(p(pNegative)))./sum(sum(pNegative));
            retPhase(:,:,iter)=p;
            %                     pause()
        end
        figure(1),
        for kk=1:6
            subplot(3,3,kk),imagesc(squeeze(retPhase(:,:,round(size(tomogMap,3)/6*kk))),[-3 3]),axis image,axis off
        end
        subplot(3,3,[7 9]),imagesc(squeeze(retPhase(:,end/2,:)),[-3 3])
        colormap('jet')
        fileName=strcat('Field_',sampleList(sampleNum).name,'_rev2');
        cd(spath)
        cd(batchList(batchNum).name)
        cd('field_retrieval');
        save(strcat(fileName,'.mat'),'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
        print('-dpng',strcat(fileName,'.png'))
    end
end
%% Tomogram Reconstruction
clear all
close all
% Angle load
spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK';
cd(spath);
batchList=dir('batch*');
for batchNum=1:length(batchList)
    cd(spath);
    cd(batchList(batchNum).name)
    cd('field_retrieval')
    sampleList=dir('Field*rev2.mat');
    for sampleNum=1:length(sampleList)
        load(sampleList((sampleNum)).name);
        [xx, yy, frame]=size(retPhase);
        crop_size=xx;
        f_dx2=f_dx-mean(f_dx(:)); % subtract maxpoint
        f_dy2=f_dy-mean(f_dy(:));
        original_size=xSize;
        excludeFrame=[];
        temp=mean(squeeze(mean(abs(retPhase),1)));
        %                                 plot(temp,'or'),hold on
        excludeFrame=[excludeFrame, find(abs(temp)>1.5)];
        temp=temp-circshift(temp,1);
        %                                 plot(temp,'og'),hold on
        excludeFrame=[excludeFrame,  find(abs(temp)>0.1)];
        %                                     end
        %                                 end
        for kkk=1:frame
            p2=squeeze(retPhase(:,:,kkk));
            if sum(isnan(p2(:)))
                excludeFrame=[excludeFrame,kkk];
            end
        end
        % figure,
        % for kk=1:size(retPhase,3)
        %     imagesc(retPhase(:,:,kk),[-3 3]),axis image,colormap('jet')
        %     title(num2str(kk))
        %     pause()
        % end
        excludeFrame=unique(sort([excludeFrame]));
        n_m=1.337;
        n_s=1.337+0.04;
        ZP=round(1.2*xx/2)*2;
        crop_factor=crop_size/original_size;
        res2=res/crop_factor;
        padd_factor=ZP/crop_size;
        kres=1/(res*ZP)*crop_factor;
        f_dx2=f_dx2*padd_factor;f_dy2=f_dy2*padd_factor;
        k0_x=kres*f_dx2; % for actual value, multiply resolution
        k0_y=kres*f_dy2;
        k0=1/lambda;
        k0_z=real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2)); % magnitude of absolute value is k0
        %%
        % Range=[54 539];
        ZP2=round(512*1.0/2)*2;
        ZP3=round(256*1.0/2)*2;
        res3=res2*ZP/ZP2;
        res4=res2*ZP/ZP3;
        frameList=(1:frame);
        frameList(excludeFrame)=[];
        IterNum=100;
        TomoParam=struct('n_m',n_m,'ZP',ZP,'ZP2',ZP2,'ZP3',ZP3,'xx',xx,'yy',yy,...
            'f_dx2',f_dx2,'f_dy2',f_dy2,'NA',NA,'lambda',lambda,...
            'k0',k0,'k0_x',k0_x,'k0_y',k0_y,'k0_z',k0_z,'kres',kres,'frameList',frameList,...
            'res2',res2,'res3',	res3,'res4',res4);
        %% Tomogram Reconstruction
        [Reconimg, ORytov]=ODTReconstruction(retAmplitude,retPhase,TomoParam);
        Reconimg=ODTIteration(gpuArray((Reconimg)),ORytov,TomoParam,IterNum);
        ZP4=round(size(Reconimg,1)/(2*padd_factor));
        ZP5=round(size(Reconimg,3)/(4*padd_factor));
        Reconimg=real(Reconimg(end/2-ZP4+1:end/2+ZP4,end/2-ZP4+1:end/2+ZP4,end/2-ZP5+1:end/2+ZP5));
        ZP4=round(size(Reconimg,1));
        ZP5=round(size(Reconimg,3));
        figure(1),
        %%
        subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-0,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
        subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-0,:)))'),[1.337-0.005 n_s]),axis image,colorbar
        %         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
        subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+1)))),[1.337-0.005 n_s]),axis image,colorbar
        subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
        colormap('jet')
        %%
        fileName=strcat('Tomogram_',sampleList((sampleNum)).name);
        save(strcat(fileName,'.mat'),'Reconimg','res3','res4','lambda','excludeFrame');
        print('-dpng',strcat(fileName,'.png'))
        clf
    end
end
%%
clear all
close all

path='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK/batch02_confluent_1day_differentiation/';
cd(path)
cd('field_retrieval')

n_m=1.337;
n_s=1.372;

tomoList2=dir(strcat('Tomogram_Field_sample',sprintf('%03d',1),'_meshScan*.mat_rev2.mat.mat'));
load(tomoList2(1).name);

tile_rows=sqrt(length(tomoList2));
tile_cols=sqrt(length(tomoList2));
% Read all tiles into a cell array
tiles = cell(tile_rows, tile_cols);
idx = 1;
for r = 1:tile_rows
    for c = 1:tile_cols
        if mod(r,2)==1 % left-to-right
            col = tile_cols - c + 1;
        else % right-to-left (bidirectional scan)
            col = c;
        end
        load(tomoList2(idx).name);
        tiles{r,col} =max(Reconimg(:,:,end/2-end/4:end/2+end/4),[],3);
        idx = idx + 1;
    end
end

% idx=1;
% figure,
% for r = 1:tile_rows
%     for c = 1:tile_cols
%         subplot(7,7,idx),imagesc(tiles{r,c},[1.34 1.36]),axis image,axis off
%         idx=idx+1;
%     end
% end

tile_size=size(Reconimg);
% Initialize stitched image canvas (rough estimate, adjusted later)
stitched = zeros(tile_rows * tile_size(1), tile_cols * tile_size(2));
positions = zeros(tile_rows, tile_cols, 2); % y, x shifts

% Set the first tile as anchor
stitched(1:tile_size(1), 1:tile_size(2)) = tiles{1,1};
positions(1,1,:) = [1, 1];

% Stitch row by row using normalized cross-correlation
for r = 1:tile_rows
    for c = 1:tile_cols
        if r == 1 && c == 1
            continue
        end
        if c > 1
            ref = (tiles{r, c-1});
            curr = (tiles{r, c});
            corrVal = xcorr2(gradient(ref), gradient(curr));
            corrVal(end/2-tile_size(1)/2:end/2+tile_size(1)/2,end/2-tile_size(2)/2:end/2+tile_size(2)/2)=0;
            [~, imax] = max(abs(corrVal(:)));
            [ypeak, xpeak] = ind2sub(size(corrVal), imax);
            corr_offset = [ypeak - size(ref,1), xpeak - size(ref,2)];
            offset = corr_offset;
            
            ref_pos = squeeze(positions(r, c-1, :))';
            new_pos = ref_pos + offset;
        else
            ref = tiles{r-1, c};
            curr = tiles{r, c};
            corrVal = xcorr2(gradient(ref), gradient(curr));
            corrVal(end/2-tile_size(1)/2:end/2+tile_size(1)/2,end/2-tile_size(2)/2:end/2+tile_size(2)/2)=0;
            [maxVal, imax] = max(abs(corrVal(:)));
            [ypeak, xpeak] = ind2sub(size(corrVal), imax);
            corr_offset = [ypeak - size(ref,1), xpeak - size(ref,2)];
            offset = corr_offset;
            
            ref_pos = squeeze(positions(r-1, c, :))';
            new_pos = ref_pos + offset;
        end
        positions(r, c, :) = new_pos;
    end
end

% Normalize positions and create stitched canvas
min_y = min(positions(:,:,1), [], 'all');
min_x = min(positions(:,:,2), [], 'all');
positions(:,:,1) = positions(:,:,1) - min_y + 1;
positions(:,:,2) = positions(:,:,2) - min_x + 1;
max_y = max(positions(:,:,1), [], 'all') + tile_size(1);
max_x = max(positions(:,:,2), [], 'all') + tile_size(2);
stitched = zeros(max_y, max_x, 'like', tiles{1,1});

% Place tiles
for r = 1:tile_rows
    for c = 1:tile_cols
        y = round(positions(r,c,1));
        x = round(positions(r,c,2));
        tile = tiles{r,c};
        stitched(y:y+tile_size(1)-1, x:x+tile_size(2)-1) = tile;
    end
end

% Display result

ZP4=round(size(stitched,1));
ZP5=round(size(stitched,3));
figure,imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,conv2(stitched,fspecial('disk',1),'same'),[n_m n_s]),axis image
xlabel('\mum');
ylabel('\mum');

save('stitched_Tomogram.mat','stitched','res3','res4');
%%
stitched(find(stitched==0))=n_m;

[H, W, ~] = size(stitched);

% Create figure and axes, using imagesc to display the image.
fig = figure('Color','black');
ax = axes('Parent', fig);
imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,conv2(stitched,fspecial('disk',1),'same'), 'Parent', ax);

axis(ax, 'image');   % Preserve aspect ratio
colormap('bone');      % Change colormap as needed
ax.CLim=[1.342 1.362];
xlabel('x (\mum)');
ylabel('y (\mum)');
c=colorbar;
c.Color='w';
c.Position=[0.902365185674867,0.115084247291850,0.035030103995621,0.804580547588513];
ax.XColor='w';
ax.YColor='w';

set(gcf,'Position',[680,580,609,518]);
%% Tomogram reconsruction version 2
clear all
close all
% Angle load

spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK/batch02_confluent_1day_differentiation';
cd(spath);
cd('field_retrieval')
sampleList=dir('Field*mesh*rev2.mat');

for sampleNum=1:length(sampleList)
    load(sampleList((sampleNum)).name);
    [~, ~, frame]=size(retPhase);
    
    f_dx2=f_dx-mean(f_dx(:)); % subtract maxpoint
    f_dy2=f_dy-mean(f_dy(:));
    original_size=xSize;
    
    excludeFrame=[];
    temp=mean(squeeze(mean(abs(retPhase),1)));
    %     plot(temp,'or'),hold on
    excludeFrame=[excludeFrame, find(abs(temp)>1.5 )];
    temp=temp-circshift(temp,1);
    %     plot(temp,'og'),hold on
    excludeFrame=[excludeFrame,  find(abs(temp)>0.1)];
    % end
    %     end
    % end
    for kkk=1:frame
        p2=squeeze(retPhase(:,:,kkk));
        if sum(isnan(p2(:)))
            excludeFrame=[excludeFrame,kkk];
        end
    end
    excludeFrame=unique(sort(excludeFrame));
    
    n_m=1.337;
    n_s=1.337+0.04;
    
    [imgH, imgW, ~] = size(retPhase);
    
    % Parameters
    numTilesX = 2;
    numTilesY = 2;
    
    % Base tile size (without overlap)
    baseW = 2*floor(imgW / numTilesX/2);
    baseH = 2*floor(imgH / numTilesY/2);
    
    % Overlap (20%)
    overlapRatio = 0.0;
    tileW = 2*round((baseW * (1 + overlapRatio))/2);  % 512 * 1.2 = 614
    tileH = 2*round((baseH * (1 + overlapRatio))/2);
    
    % Stride stays at base tile size
    strideX = baseW;
    strideY = baseH;
    
    % Initialize tile storage
    phaseTiles = {};
    amplitudeTiles = {};
    tileNum = 1;
    
    for row = 0:(numTilesY - 1)
        for col = 0:(numTilesX - 1)
            % Top-left corner
            x_start = col * strideX + 1;
            y_start = row * strideY + 1;
            
            % Make sure we don't exceed image boundary
            x_end = min(x_start + tileW - 1, imgW);
            y_end = min(y_start + tileH - 1, imgH);
            
            x_start = max(1, x_end - tileW + 1);
            y_start = max(1, y_end - tileH + 1);
            
            % Extract tile
            phaseTile = retPhase(y_start:y_end, x_start:x_end, :);
            amplitudeTile = retAmplitude(y_start:y_end, x_start:x_end, :);
            phaseTiles{tileNum} = phaseTile;
            amplitudeTiles{tileNum} = amplitudeTile;
            
            % Optional: save
            % imwrite(tile, sprintf('tile_%02d.png', tileNum));
            
            tileNum = tileNum + 1;
        end
    end
    %     figure
    %     for mm=1:length(phaseTiles)
    %         retPhase2=phaseTiles{mm};
    %
    %         subplot(4,4,mm),imagesc(retPhase2(:,:,49),[-3 3]),axis image,colormap('jet'),axis off
    %     end
    ReconimgTiles = {};
    tileNum = 1;
    for mm=1:length(phaseTiles)
        retPhase2=phaseTiles{mm};
        retAmplitude2=amplitudeTiles{mm};
        
        [xx, yy, frame]=size(retPhase);
        [xx2, yy2, frame]=size(retPhase2);
        
        crop_size=imgH;
        ZP=round(1.2*xx2/2)*2;
        crop_factor=crop_size/original_size;
        res2=res/crop_factor;
        padd_factor=ZP/xx;
        padd_factor2=ZP/xx2;
        kres=1/(res*ZP)*crop_factor;
        f_dx3=f_dx2*padd_factor;f_dy3=f_dy2*padd_factor;
        k0_x=kres*f_dx3; % for actual value, multiply resolution
        k0_y=kres*f_dy3;
        k0=1/lambda;
        k0_z=real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2)); % magnitude of absolute value is k0
        %%
        % Range=[54 539];
        ZP2=round(512*1.0/2)*2;
        ZP3=round(256*1.0/2)*2;
        res3=res2*ZP/ZP2;
        res4=res2*ZP/ZP3;
        frameList=(1:frame);
        frameList(excludeFrame)=[];
        IterNum=10;
        TomoParam=struct('n_m',n_m,'ZP',ZP,'ZP2',ZP2,'ZP3',ZP3,'xx',xx2,'yy',yy2,...
            'f_dx2',f_dx3,'f_dy2',f_dy3,'NA',NA,'lambda',lambda,...
            'k0',k0,'k0_x',k0_x,'k0_y',k0_y,'k0_z',k0_z,'kres',kres,'frameList',frameList,...
            'res2',res2,'res3',	res3,'res4',res4);
        %% Tomogram Reconstruction
        [Reconimg, ORytov]=ODTReconstruction_scale(retAmplitude2,retPhase2,TomoParam);
        %         Reconimg=(gather(fftshift(Reconimg,3)));
        Reconimg=ODTIteration(gpuArray((Reconimg)),ORytov,TomoParam,IterNum);
        %         Reconimg=ODTIteration_scale(gpuArray((Reconimg)),ORytov,TomoParam,IterNum);
        
        ZP4=round(size(Reconimg,1)/(2*padd_factor2));
        ZP5=round(size(Reconimg,3)/(4*padd_factor2));
        Reconimg=real(Reconimg(end/2-ZP4+1:end/2+ZP4,end/2-ZP4+1:end/2+ZP4,end/2-ZP5+1:end/2+ZP5));
        ReconimgTiles{mm}=Reconimg;
    end
    
    %         figure
    %         for mm=1:length(phaseTiles)
    %             retPhase2=ReconimgTiles{mm};
    %
    %             subplot(numTilesX,numTilesY,mm),imagesc(max(retPhase2,[],3),[1.34 n_s]),axis image,colormap('jet'),axis off
    %         end
    
    tileH = size(ReconimgTiles{1}, 1);
    tileW = size(ReconimgTiles{1}, 2);
    numChannels = size(ReconimgTiles{1}, 3);  % For grayscale use 1, RGB use 3
    
    % Initialize the stitched image
    stitchedImage = zeros(numTilesX * tileH, numTilesY * tileW, numChannels, 'like', ReconimgTiles{1});
    
    % Fill in the stitched image
    tileIdx = 1;
    for row = 0:(numTilesX - 1)
        for col = 0:(numTilesY - 1)
            yStart = row * tileH + 1;
            yEnd = (row + 1) * tileH;
            xStart = col * tileW + 1;
            xEnd = (col + 1) * tileW;
            
            stitchedImage(yStart:yEnd, xStart:xEnd, :) = ReconimgTiles{tileIdx};
            tileIdx = tileIdx + 1;
        end
    end
    
    ZP4=round(size(stitchedImage,1));
    ZP5=round(size(stitchedImage,3));
    figure(1),
    %%
    subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(stitchedImage(end/2+-0,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
    subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(stitchedImage(:,end/2+-0,:)))'),[1.337-0.005 n_s]),axis image,colorbar
    %         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
    subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(stitchedImage(:,:,end/2+1)))),[1.337-0.005 n_s]),axis image,colorbar
    subplot(224),imagesc(max(real(stitchedImage(:,:,end/2-end/4:end/2+end/4)),[],3),[1.337-0.005 n_s]),axis image, axis off
    colormap('jet')
    
    %%
    fileName=strcat('Tomogram_',sampleList((sampleNum)).name,'_v3');
    save(strcat(fileName,'.mat'),'stitchedImage','res3','res4','lambda');
    print('-dpng',strcat(fileName,'.png'))
    clf
end

%%
clear all
close all

path='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK/batch02_confluent_1day_differentiation/';
cd(path)
cd('field_retrieval')

n_m=1.337;
n_s=1.372;

tomoList2=dir(strcat('Tomogram_Field_sample',sprintf('%03d',1),'_meshScan*.mat_v3.mat'));
load(tomoList2(1).name);

tile_rows=sqrt(length(tomoList2));
tile_cols=sqrt(length(tomoList2));
% Read all tiles into a cell array
tiles = cell(tile_rows, tile_cols);
idx = 1;
for r = 1:tile_rows
    for c = 1:tile_cols
        if mod(r,2)==1 % left-to-right
            col = tile_cols - c + 1;
        else % right-to-left (bidirectional scan)
            col = c;
        end
        load(tomoList2(idx).name);
        tiles{r,col} =max(stitchedImage(:,:,end/2-end/4:end/2+end/4),[],3);
        idx = idx + 1;
    end
end

% idx=1;
% figure,
% for r = 1:tile_rows
%     for c = 1:tile_cols
%         subplot(7,7,idx),imagesc(tiles{r,c},[1.34 1.36]),axis image,axis off
%         idx=idx+1;
%     end
% end

tile_size=size(stitchedImage);
% Initialize stitched image canvas (rough estimate, adjusted later)
stitched = zeros(tile_rows * tile_size(1), tile_cols * tile_size(2));
positions = zeros(tile_rows, tile_cols, 2); % y, x shifts

% Set the first tile as anchor
stitched(1:tile_size(1), 1:tile_size(2)) = tiles{1,1};
positions(1,1,:) = [1, 1];

% Stitch row by row using normalized cross-correlation
for r = 1:tile_rows
    for c = 1:tile_cols
        if r == 1 && c == 1
            continue
        end
        if c > 1
            ref = (tiles{r, c-1});
            curr = (tiles{r, c});
            corrVal = xcorr2(gradient(ref), gradient(curr));
            corrVal(end/2-tile_size(1)/2:end/2+tile_size(1)/2,end/2-tile_size(2)/2:end/2+tile_size(2)/2)=0;
            [~, imax] = max(abs(corrVal(:)));
            [ypeak, xpeak] = ind2sub(size(corrVal), imax);
            corr_offset = [ypeak - size(ref,1), xpeak - size(ref,2)];
            offset = corr_offset;
            
            ref_pos = squeeze(positions(r, c-1, :))';
            new_pos = ref_pos + offset;
        else
            ref = tiles{r-1, c};
            curr = tiles{r, c};
            corrVal = xcorr2(gradient(ref), gradient(curr));
            corrVal(end/2-tile_size(1)/2:end/2+tile_size(1)/2,end/2-tile_size(2)/2:end/2+tile_size(2)/2)=0;
            [maxVal, imax] = max(abs(corrVal(:)));
            [ypeak, xpeak] = ind2sub(size(corrVal), imax);
            corr_offset = [ypeak - size(ref,1), xpeak - size(ref,2)];
            offset = corr_offset;
            
            ref_pos = squeeze(positions(r-1, c, :))';
            new_pos = ref_pos + offset;
        end
        positions(r, c, :) = new_pos;
    end
end

% Normalize positions and create stitched canvas
min_y = min(positions(:,:,1), [], 'all');
min_x = min(positions(:,:,2), [], 'all');
positions(:,:,1) = positions(:,:,1) - min_y + 1;
positions(:,:,2) = positions(:,:,2) - min_x + 1;
max_y = max(positions(:,:,1), [], 'all') + tile_size(1);
max_x = max(positions(:,:,2), [], 'all') + tile_size(2);
stitched = zeros(max_y, max_x, 'like', tiles{1,1});

% Place tiles
for r = 1:tile_rows
    for c = 1:tile_cols
        y = round(positions(r,c,1));
        x = round(positions(r,c,2));
        tile = tiles{r,c};
        stitched(y:y+tile_size(1)-1, x:x+tile_size(2)-1) = tile;
    end
end

% Display result

ZP4=round(size(stitched,1));
ZP5=round(size(stitched,3));
figure,imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,conv2(stitched,fspecial('disk',1),'same'),[n_m n_s]),axis image
xlabel('\mum');
ylabel('\mum');

save('stitched_Tomogram_v3.mat','stitched','res3','res4');

%% Visualization
clear all
close all

spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK';
cd(spath);
batchList=dir('batch*');
for batchNum=1:length(batchList)
    cd(spath);
    cd(batchList(batchNum).name)
    cd('field_retrieval')
    sampleList=dir('Tomogram_Field_sample*_Tomog.mat_rev2.mat.mat');
    figure,
    for sampleNum=1:length(sampleList)
        load(sampleList((sampleNum)).name);
        
        subplot(3,4,sampleNum)
        imagesc(max(Reconimg,[],3),[1.337 1.37]),axis image,colormap('jet'),axis off
    end
end
%%
n_m=1.337;
n_s=1.37;

spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK';
cd(spath);
batchList=dir('batch*');
for batchNum=1:length(batchList)
    cd(spath);
    cd(batchList(batchNum).name)
    cd('field_retrieval')
    sampleList=dir('Tomogram_Field_sample*_Tomog.mat_rev2.mat.mat');
    for sampleNum=1:length(sampleList)
        load(sampleList((sampleNum)).name);
        
        ZP4=round(size(Reconimg,1));
        ZP5=round(size(Reconimg,3));
        figure,
        %%
        subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-0,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
        subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-0,:)))'),[1.337-0.005 n_s]),axis image,colorbar
        %         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
        subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+1)))),[1.337-0.005 n_s]),axis image,colorbar
        subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
        colormap('jet')
    end
end

%%
spath='/beegfs/home/ralajan/matlab/20260203_Cecile_MDCK';
cd(spath);
batchList=dir('batch*');
batchNum=2;
cd(spath);
cd(batchList(batchNum).name)
cd('field_retrieval')
sampleList=dir('Tomogram_Field_sample*_Tomog.mat_rev2.mat.mat');
sampleNum=6;
load(sampleList((sampleNum)).name);

ZP4=round(size(Reconimg,1));
ZP5=round(size(Reconimg,3));
figure,
%%
subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-19,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-5,:)))'),[1.337-0.005 n_s]),axis image,colorbar
%         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+3)))),[1.337-0.005 n_s]),axis image,colorbar
subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
colormap('jet')

%%
sampleNum=8;
load(sampleList((sampleNum)).name);

ZP4=round(size(Reconimg,1));
ZP5=round(size(Reconimg,3));
figure,
%%
subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-20,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-4,:)))'),[1.337-0.005 n_s]),axis image,colorbar
%         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+3)))),[1.337-0.005 n_s]),axis image,colorbar
subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
colormap('jet')

%%
sampleNum=10;
load(sampleList((sampleNum)).name);

ZP4=round(size(Reconimg,1));
ZP5=round(size(Reconimg,3));
figure,
%%
subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-23,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-4,:)))'),[1.337-0.005 n_s]),axis image,colorbar
%         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+8)))),[1.337-0.005 n_s]),axis image,colorbar
subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
colormap('jet')
%%
sampleNum=11;
load(sampleList((sampleNum)).name);

ZP4=round(size(Reconimg,1));
ZP5=round(size(Reconimg,3));
figure,
%%
subplot(221),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(end/2+-3,:,:)))'),[1.337-0.005 n_s]),axis image,colorbar
subplot(222),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP5)-ZP5/2)*res4,(real(squeeze(Reconimg(:,end/2+-4,:)))'),[1.337-0.005 n_s]),axis image,colorbar
%         subplot(4,3,12),imagesc(real(squeeze(Reconimg(:,:,end/2+z))),[n_m-0.01 1.47]),axis image,colorbar,axis off
subplot(223),imagesc(((1:ZP4)-ZP4/2)*res3,((1:ZP4)-ZP4/2)*res3,(real(squeeze(Reconimg(:,:,end/2+7)))),[1.337-0.005 n_s]),axis image,colorbar
subplot(224),imagesc(max(real(Reconimg),[],3),[1.337-0.005 n_s]),axis image, axis off
colormap('jet')
