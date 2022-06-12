clear all;
clf
load U_3D.mat;
file_out='u_z.grd';
xb=[-118 -117]; % longtitude window : Ridgecrest
yb=[35.2 36];   % latitude window

%return

%grdwrite2(x(1,:),y(:,1),ux,'ux.grd');
%grdwrite2(x(1,:),y(:,1),uy,'uy.grd');
grdwrite2(x(1,:),y(:,1),uz,file_out);

