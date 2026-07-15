%% Angle load
% spath is provided by main.sh at invocation
bglist = dir(fullfile(spath, 'bg*_Tomog.mat'));
f_dx = [];
f_dy = [];
load(fullfile(spath, bglist(1).name));
img = double(squeeze(tomogMap(:,:,round(size(tomogMap,3)/3)-1)));
img = squeeze(img);
ii = length(img);
ZP = ii;
r = round(ZP*res*NA/lambda)+20;
yr = r;
for bgNum = 1:size(tomogMap,3)
    img = double(squeeze(tomogMap(:,:,bgNum)));
    [xSize ySize] = size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize);
    [f_x,f_y] = find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49)))));
    f_dx = [f_dx; f_x];
    f_dy = [f_dy; f_y];
end
mi = mean(f_dx);
mj = mean(f_dy);
mi = round(mi-ii/2-1); mj = round(mj-ii/2-1);
c1mask = ~(mk_ellipse(r-20,yr-20,ZP,ZP));
c3mask = circshift(c1mask,[mi mj]);
Fbg = zeros(round(2*r),round(2*r),size(tomogMap,3),'single');
for bgNum = 1:size(tomogMap,3)
    img = double(squeeze(tomogMap(:,:,bgNum)));
    [xSize ySize] = size(img);
    Fimg = fftshift(fft2(img))/(xSize*ySize);
    Fimg = Fimg.*c3mask;
    Fimg = circshift(Fimg,-[round(mi) round(mj)]);
    Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
    sizeFimg = length(Fimg);
    Fimg = ifft2(ifftshift(Fimg))*(sizeFimg^2);
    Fbg(:,:,bgNum) = Fimg;
end

%%
outDir = fullfile(spath, 'field_retrieval');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

sampleList = dir(fullfile(spath, 'sample*_Tomog.mat'));
for sampleNum = 1:length(sampleList)
    load(fullfile(spath, sampleList(sampleNum).name));
    retPhase = zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
    retAmplitude = zeros(round(2*r)-4,round(2*r)-4,size(tomogMap,3),'single');
    f_dx = [];
    f_dy = [];
    for iter = 1:size(tomogMap,3)
        img = (double(squeeze(tomogMap(:,:,iter))));
        Fimg = fftshift(fft2(img))/(ii^2);
        [f_x,f_y] = find(Fimg==max(max(Fimg(:, round(ii*0.01):round(ii*0.49)))));
        f_dx = [f_dx; f_x];
        f_dy = [f_dy; f_y];
        Fimg = Fimg.*c3mask;
        Fimg = circshift(Fimg,-[round(mi) round(mj)]);
        Fimg = Fimg(ii/2-r+1:ii/2+r,ii/2-r+1:ii/2+r);
        sizeFimg = length(Fimg);
        Fimg = (ifft2(ifftshift(Fimg))*(sizeFimg^2));
        Fimg = Fimg./squeeze(Fbg(:,:,iter));
        %%
        FFimg = (fftshift(fft2(Fimg)));
        [tempX,tempY] = find(abs(FFimg)==max(max(abs(FFimg))));
        tempX = tempX-size(FFimg,1)/2;
        tempY = tempY-size(FFimg,2)/2;
        Fimg = ifft2(ifftshift(circshift(FFimg,-[tempX-1 tempY-1])));
        %%
        Fimg = Fimg(3:end-2,3:end-2);
        retAmplitude(:,:,iter) = abs(Fimg);
        p = (PhiShift(unwrap2(double(angle(Fimg)))));
        [p,coeffVal,~] = phaseCompensation(p,1);
        pp = (p);
        tempThresh = p-imtophat(pp,strel('disk',150));
        tempThresh = mean(tempThresh(:))+1.0;
        p2Mask = (pp>tempThresh(1));
        p2Mask = bwareaopen(p2Mask,100);
        p2Mask = ~imdilate(p2Mask,strel('disk',5));
        [p,coefficients] = phaseCompensation(p,2,p2Mask);
        pNegative = (p<0);
        p = p-sum(sum(p(pNegative)))./sum(sum(pNegative));
        retPhase(:,:,iter) = p;
    end
    figure(1),
    for kk = 1:6
        subplot(3,3,kk),imagesc(squeeze(retPhase(:,:,round(size(tomogMap,3)/6*kk))),[-3 3]),axis image,axis off
    end
    subplot(3,3,[7 9]),imagesc(squeeze(retPhase(:,end/2,:)),[-3 3])
    colormap('jet')
    fileName = strcat('Field_',sampleList(sampleNum).name);
    save(fullfile(outDir, strcat(fileName,'.mat')),'retAmplitude','retPhase','xSize','f_dx','f_dy','NA','lambda','res','ZP');
    print('-dpng', fullfile(outDir, strcat(fileName,'.png')))
end

%% Field inspection plot
fieldList = dir(fullfile(outDir, 'Field*.mat'));
figure(2), hold on
for sampleNum = 1:length(fieldList)
    load(fullfile(outDir, fieldList(sampleNum).name));
    f_dx2 = f_dx-mean(f_dx(:));
    f_dy2 = f_dy-mean(f_dy(:));
    temp = mean(squeeze(mean(abs(retPhase),1)));
    plot(temp,'or')
    temp = temp-circshift(temp,1);
    plot(temp,'og')
end
ylim([-3 3]);
print('-dpng', fullfile(outDir, 'Field_inspection.png'))
