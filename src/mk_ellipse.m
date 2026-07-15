function H = mk_ellipse(XR,YR,X,Y)
 
[XX YY]=meshgrid(1:X,1:Y);
H=(((XX-X/2)./XR).^2+((YY-Y/2)./YR).^2)>1.0;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Out = PhaseShift (Uimg,px,py,isize,jsize)
is=10;
ie=isize-10;
js=10;
je=jsize-10;
if px==0 && py==0
  px=sum((Uimg(is:ie,js)-Uimg(is:ie,je)))/((ie-is+1)*(je-js+1));
  py=sum((Uimg(is,js:je)-Uimg(ie,js:je)))/((ie-is+1)*(je-js+1));
end
[XX YY]=meshgrid(1:jsize,1:isize);
Out=Uimg+XX.*px+YY.*py;
daout = Out;
subplot(2,2,1); plot(daout(200,:)), axis ([0 640 -6 6]),title('X-axis');
subplot(2,2,2); plot(daout(:,500)), axis([0 480 -6 6]),title('Y-axis');
%status = fclose('all');
return;