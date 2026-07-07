function Out = PhiShift(Uimg)
%Uimg=Timg(1:200,1:240);
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


% figure(100),imagesc(Out)
return;