 clear all;
clf
load U_3D.mat;
file_out='u_disp.grd';
uh=sqrt(ux.^2+uy.^2);
xb=[-118 -117]; % longtitude window : Ridgecrest
yb=[35.2 36];   % latitude window
 [m,n]=size(x);
 fct=5.7e1;
 dhx=round(n/fct);
 dhy=round(m/fct);
 xsub = x(1:dhy:m,1:dhx:n);
 ysub = y(1:dhy:m,1:dhx:n);
 uysub = uy(1:dhy:m,1:dhx:n);
 uxsub = ux(1:dhy:m,1:dhx:n);
figure(1)
clf
 fs=200;
 quiver(xsub,ysub,uxsub,uysub,1e0,'k'), hold on
% set(gca,'xlim',xb,'ylim',yb), hold on
fid = fopen('cos_disp.ll','w+');
for i=1:length(xsub(:))
  if ~isnan(uxsub(i))
    len=sqrt(uxsub(i)^2+uysub(i)^2);
    azi=atan2(uxsub(i),uysub(i))/pi*180;
%   fprintf(fid,'%f %f %f %f \n',xsub(i),ysub(i),uxsub(i)/fs,uysub(i)/fs);
   fprintf(fid,'%f %f %f %f \n',xsub(i),ysub(i),azi,len/fs);
  end
end
status = fclose(fid);

%return

grdwrite2(x(1,:),y(:,1),ux,'ux.grd');
grdwrite2(x(1,:),y(:,1),uy,'uy.grd');
%grdwrite2(x(1,:),y(:,1),uh,file_out);

