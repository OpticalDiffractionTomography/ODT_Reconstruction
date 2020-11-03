function [goodimg,mdx,mdy] = PhiShiftJH(varargin)
% Original code from Dr. J. Jeong (jeongjh513@kaist.ac.kr)

% PhiShiftJH(img): Phase gradient is calculated based on the phase values at the edges of img
% img should be larger than 12-by-12 for proper operation.
% If length of any axis of img is less than 5, the error will occur.
% PhiShiftJH(img,mask): Phase gradient is calculated based on the phase values in binary mask
switch nargin
    case 1
        Uimg=varargin{1};
        [ie, je] = size(Uimg);
        b1=2;
        bsize=4;

        left=mean2(Uimg(:,b1:b1+bsize));
        right=mean2(Uimg(:, je-bsize-b1+1:je-b1+1));
        top=mean2(Uimg(b1:b1+bsize, 1:end));
        bottom=mean2(Uimg(ie-bsize-b1+1:ie-b1+1,:));

        mdx=(right-left)/je;
        mdy=(bottom-top)/ie;

        [XX YY]=meshgrid(1:je,1:ie);
        goodimg=Uimg-(XX.*mdx+YY.*mdy);
        ix=[b1+bsize, b1+bsize, je-b1-bsize, je-b1-bsize];
        iy=[b1+bsize, ie-b1-bsize, ie-b1-bsize, b1+bsize];
        boundary_mask=~poly2mask(ix,iy,ie,je);
        boundary_mask(1,:)=0;
        boundary_mask(end,:)=0;
        boundary_mask(:,1)=0;
        boundary_mask(:,end)=0;
        goodimg=goodimg-sum(sum(goodimg.*boundary_mask))./sum(sum(boundary_mask));    

    case 2;
        Uimg=varargin{1};
        BGmask=varargin{2};
        
        [yy,xx] = size(Uimg);
        
        BGmasky=BGmask(1:end-1,:);
        BGmaskx=BGmask(:,1:end-1);
        
        dBGy=diff(BGmask,1,1);
        dBGx=diff(BGmask,1,2);
        BGmasky(dBGy==-1)=0;
        BGmaskx(dBGx==-1)=0;
        
        dy_Uimg = diff(Uimg,1,1).*BGmasky;
        dx_Uimg = diff(Uimg,1,2).*BGmaskx;

        mdy = sum(sum(dy_Uimg))./sum(sum(BGmasky)); 
        mdx = sum(sum(dx_Uimg))./sum(sum(BGmaskx)); 

        [XX,YY] = meshgrid(1:xx,1:yy);
        cphi = mdy*YY+mdx*XX;
        goodimg = (Uimg - cphi);
        goodimg = goodimg-sum(sum(goodimg.*BGmask))/sum(sum(BGmask));
end