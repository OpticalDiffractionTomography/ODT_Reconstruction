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