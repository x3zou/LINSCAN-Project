clear all;
close all;

bgrey = [1 1 1];
wid=2;
fault_tr = 1; % flag for plotting the fault trace
%load Hauksson_ridgecrest_17_July.reloc;
%x=Hauksson_ridgecrest_17_July;
%X=x(:,3);
%Y=x(:,2);
%Z=x(:,4);
%plot(X, Y,'.k'), hold on

load eqn3.ross;
x=eqn3;

f=figure;
%f = figure('visible','off');
width = 800;
height = 800;
set(f,'Position',[15 15 width height])

%set(gcf,'Renderer','zbuffer') 
%set(gca,'Color','w')

X=x(:,1);
Y=x(:,2);
Z=x(:,3);

orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
x=X; y=Y;

[X1,Y1]=utm2ll(x,y,0,1);
xU=(X1-xo)*1e-3;
yU=(Y1-yo)*1e-3;
 
%plot(xU, yU,'.k'), hold on
%return

plot3(xU, yU, Z,'k.','MarkerSize',5), hold on

%load Ridgecrest1.trace;
%x=Ridgecrest1;
%X=x(:,1);
%Y=x(:,2);
%Z=0*x(:,1);
%%event = find(X> -117.9 & X < -117.2 & Y < 35.85 & Y > 35.75);
%event = find(X> -118.2 & X < -117 & Y < 36 & Y > 35);
%x=X(event);
%y=Y(event);

%[X1,Y1]=utm2ll(x,y,0,1);
%xU=(X1-xo)*1e-3;
%yU=(Y1-yo)*1e-3;

minlon=min(xU); maxlon=max(xU);
minlat=min(yU); maxlat=max(yU);

%set(gca,'xlim',[-20 -3],'ylim',[27 39],'zlim',[-11 0]);

if fault_tr==1
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');
 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');
end

load ../figs/MODEL3_Z.mat;

%NS=30;

npl=0;                 % plot slip distribution
SS=2;                 % slip scale 
ut1=[];
ut2=[];
xs=[];xv=[];yv=[];zv=[];u1=[];u2=[];u3=[];
fact=1;
Mo=0;
Moss=0;
MM=0;
xs=[];
NT=sum(Mode);
u=US(1:sum(tSm)*NT);
[U]=split_U(u,Mode);
if Mode(3)==0
 cmax=max(-U(1,:));
 cmin=min(-U(1,:));
else
 cmax=max(U(3,:));
 cmin=min(U(3,:));
% cmin=0;
end 

pos=[];
pos(1,:)=[0.1 0.3 0 0];
h=[];
avc=0; av_sl=[];
flen(1)=0;
M_T=zeros(1,6);
for i= 1:NS
    XB=[];
    YB=[];
    XO=[];
    YO=[];
 if i > 1, MM=MI(i-1); end
  ias=(1:MI(i))+MM*(i-1);
  if src_type(i) == 'O' | src_type(i) == 'o'
    asO=as(ias);
    npl=npl+1;
    [depth,xc,yc,XB,YB,XO,YO]=grid_plane_zeyu(asO,Lp,Wp,XB,YB,XO,YO,dL,dW,nL,npl);
    k=sum(tSm(1:npl))*NT;
    u=US(k+1:k+sum(tSm(npl+1:npl+1))*NT);
    [U]=split_U(u,Mode)/fact;
    if Mode(3)==0, U(3,:)=U(1,:); end
    if sum(tSm) > 4
     xr=[];
     yr=[];
     xmi=Inf;
     xma=-Inf;
     ymi=Inf;
     yma=-Inf;
     av_sl=[av_sl zeros(1,Wp(npl))];
     k=0;
     for jc=1:Wp(npl)
      Yc=-depth(1,jc)/abs(sin(asO(6)));
      for ic=1:nL(npl,jc)
       k=k+1;
       xb=XB((k-1)*5+1:(k-1)*5+4);
       yb=YB((k-1)*5+1:(k-1)*5+4);
       dt=-depth(ic,jc)-dW(npl,jc)*sin(asO(6));
       zb=[dt -depth(ic,jc) -depth(ic,jc) dt];
       xv=[xv 0.5*(xb(1)+xb(3))];
       yv=[yv 0.5*(yb(1)+yb(3))];
       zv=[zv 0.5*(dt-depth(ic,jc))];
       x1=U(2,k)*cos(asO(6)); y1=U(1,k);
       x2=x1*cos(asO(7))+y1*sin(asO(7));
       y2=-x1*sin(asO(7))+y1*cos(asO(7));
       u1=[u1 x2];
       u2=[u2 y2];
       u3=[u3 U(2,k)*sin(asO(6))];
       ut=sqrt(U(1,k)^2+U(2,k)^2);
       if U(1,k) > 0, U(1,k)=-U(1,k); end;
%if mean(yb) > 28 & mean(yb) < 39 & mean(xb) < 0

        if Mode(3) == 0
         h=fill3(xb,yb,zb,-U(1,k),'EdgeColor','k','LineWidth',1,'FaceAlpha',0.5,'EdgeAlpha',1); hold on  %
%         h=fill3(xb,yb,zb,-U(1,k),'EdgeColor','k','LineWidth',4,'EdgeAlpha',1); hold on  %
        else
         fill3(xb,yb,zb,U(3,k),'EdgeColor','k','FaceAlpha',0.5); hold on  %
        end
%       plot3(xc,yc,-depth(ic,jc),'ko'), hold on
%end
%set(h,'edgecolor','k','LineWidth',4);
      end
     end
%     clf
     avc=avc+Wp(npl);

     flen(i+1)=flen(i)+xma-xmi;
    end
  end
end

   caxis([0 350]);

i1=npl+1;
%return

% hypocenter 6.4
x1=-117.506; y1=35.705; z1=-10.5;
[X1,Y1]=utm2ll(x1,y1,0,1);
x1=(X1-xo)*1e-3;
y1=(Y1-yo)*1e-3;
% hypocenter 7.1
x1=-117.599; y1=35.770; z1=-8;
[X1,Y1]=utm2ll(x1,y1,0,1);
x1=(X1-xo)*1e-3;
y1=(Y1-yo)*1e-3;
[x,y,z] = sphere(100); hold on
R=1; % sphere radius, km
hs1 = surf(x*R+x1,y*R+y1,z*R+z1,'FaceAlpha',0.4); 
q1 = get(hs1);
set(hs1, 'FaceColor', [0 1 1]);

%return

% hypocenter 7.1
x1=-117.599; y1=35.770; z1=-8;
[X1,Y1]=utm2ll(x1,y1,0,1);
x1=(X1-xo)*1e-3;
y1=(Y1-yo)*1e-3;
%[x,y,z] = sphere(100); hold on
R=1; % sphere radius, km
%hs2 = surf(x*R+x1,y*R+y1,z*R+z1,'FaceAlpha',0.4,'edgecolor','k','EdgeAlpha',1); shading flat
%q2 = get(hs2);
%set(hs2, 'FaceColor', [1 0 0]);
text(-8,40,0, 'North', 'FontSize',20,'Rotation',+67)
grid on;

colormap jet;
set(gcf,'color','w');
set(gca,'Projection','Perspective','color','w');

axis('equal'), hold on
set(gca,'xlim',[-30 maxlon],'ylim',[minlat maxlat],'zlim',[-25 0]);
axis vis3d, hold on
axis manual, hold on
 view(-67,39)
 print('-djpeg90','-r600','mod_3d')

return

videoFilename = fullfile(pwd,'slip_movie.mp4');
%videoFilename = fullfile(pwd,'slip_movie.avi');
%if exist(videoFilename,'file')
%    delete(videoFilename)
%endset(gca,'Color','k')

vidfile = VideoWriter(videoFilename,'MPEG-4');
%vidfile = VideoWriter(videoFilename);
vidfile.Quality = 100;
%vidfile.Duration = 30;
vidfile.FrameRate = 15;

%if exist(vidfile,'file')
%    delete(vidfile)
%end

open(vidfile);

j=36;
%for i=0:1:360
%for i=-80:280
for i=-87:-87
%  if i < 130 
%   j=j+0.2;
%  else
%   j=j-0.2;
%  end
 view(i,j)
 drawnow
 print('-djpeg90','-r600','mod_3d')

 frame = getframe(gcf);
 writeVideo(vidfile, frame);
end


close(vidfile)


