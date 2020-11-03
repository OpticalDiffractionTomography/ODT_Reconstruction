function Out = PhiShift(Uimg)
% Phase gradient is calculated based on the phase values at the edges of Uimg
[ie, je] = size(Uimg);
b1=2;
bsize=4;

left=mean2(Uimg(:,b1:b1+bsize));
right=mean2(Uimg(:, je-bsize-b1+1:je-b1+1));
top=mean2(Uimg(b1:b1+bsize, 1:end));
bottom=mean2(Uimg(ie-bsize-b1+1:ie-b1+1,:));

px=(right-left)/je;
py=(bottom-top)/ie;

[XX YY]=meshgrid(1:je,1:ie);
Out=Uimg-(XX.*px+YY.*py);
return;