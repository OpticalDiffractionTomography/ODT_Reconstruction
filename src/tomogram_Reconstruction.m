%% Tomogram Reconstruction
% spath and addpath are provided by main.sh at invocation
cd(spath);
batchList=dir('batch*');
for batchNum=1:length(batchList)
    cd(spath);
    cd(batchList(batchNum).name)
    cd('field_retrieval')
    sampleList=dir('Field*.mat');
    for sampleNum=1:length(sampleList)
        load(sampleList((sampleNum)).name);
        [xx, yy, frame]=size(retPhase);
        crop_size=xx;
        f_dx2=f_dx-mean(f_dx(:));
        f_dy2=f_dy-mean(f_dy(:));
        original_size=xSize;
        excludeFrame=[];
        temp=mean(squeeze(mean(abs(retPhase),1)));
        excludeFrame=[excludeFrame, find(abs(temp)>1.5)];% variable input
        temp=temp-circshift(temp,1);
        excludeFrame=[excludeFrame,  find(abs(temp)>0.1)];% variable input

        for kkk=1:frame
            p2=squeeze(retPhase(:,:,kkk));
            if sum(isnan(p2(:)))
                excludeFrame=[excludeFrame,kkk];
            end
        end

        excludeFrame=unique(sort([excludeFrame]));
        n_m=1.337;
        n_s=1.337+0.04;
        ZP=round(1.2*xx/2)*2;
        crop_factor=crop_size/original_size;
        res2=res/crop_factor;
        padd_factor=ZP/crop_size;
        kres=1/(res*ZP)*crop_factor;
        f_dx2=f_dx2*padd_factor;f_dy2=f_dy2*padd_factor;
        k0_x=kres*f_dx2; 
        k0_y=kres*f_dy2;
        k0=1/lambda;
        k0_z=real(sqrt((n_m*k0)^2-(k0_x).^2-(k0_y).^2)); 
        %%
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