clear all;
close all;

fault_tr = 1; % flag for plotting the fault trace
bgrey = [1 1 1];

%load ../figs/MODEL3_Z.mat;
load MODEL3_Z.mat;

fid1 = fopen('strike_Mo.dat','wt');
fmt = ' %10.7f %10.7f \n';
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
       strike=asO(7)/pi*180;
       u1=[u1 x2];
       u2=[u2 y2];
       u3=[u3 U(2,k)*sin(asO(6))];
       if U(1,k) > 0, U(1,k)=-U(1,k); end;
        if Mode(3) == 0
         fill3(xb,yb,zb,-U(1,k),'EdgeColor','k'); hold on  %
        else
         fill3(xb,yb,zb,U(3,k),'EdgeColor','k'); hold on  %
        end
%       plot3(xc,yc,-depth(ic,jc),'ko'), hold on
       Mo1=sqrt(U(1,k)^2+U(2,k)^2)*dL(npl,jc)*dW(npl,jc);
       Mo=Mo+sqrt(U(1,k)^2+U(2,k)^2)*dL(npl,jc)*dW(npl,jc);
       Moss=Moss+U(1,k)*dL(npl,jc)*dW(npl,jc);
       av_sl(avc+jc)=av_sl(avc+jc)+abs(U(1,k));
       fprintf(fid1,fmt, [strike, Mo1]);
      end
     end
%     clf
     avc=avc+Wp(npl);

     flen(i+1)=flen(i)+xma-xmi;
    end
  end
end
fclose(fid1);

i1=npl+1;
%return
axis('equal'), hold on
%set(gca,'xlim',[-20 -3],'ylim',[27 39],'zlim',[-11 0]);
set(gca,'xlim',[-30 10],'ylim',[-10 65],'zlim',[-25 0]);

fprintf('\n Total seismic moment= %e \n',Mo*3.3e14);
fprintf('\n Total moment magnitude= %e \n',(log10(Mo*3.3e14)-9.1)/1.5);
fprintf('\n strike-slip seismic moment= %e \n',Moss*3.3e14);

return

h(i1)=subplot(1,1,1);
%axis manual
set(gca,'Projection','Perspective','color',bgrey);
%set(gca,'XLim',[-15 20],'YLim',[-25 20],'ZLim',[-20 0]);
%cmin=1e9;
%cmax=-1e9;
npl=0;                 % plot slip distribution
SS=200;                 % slip scale 
ut1=[];
ut2=[];
MM=0;
xs=[];xv=[];yv=[];zv=[];u1=[];u2=[];u3=[];

for i= 1:NS
 if i > 1, MM=MI(i-1); end
  ias=(1:MI(i))+MM*(i-1);
    asO=as(ias);
    npl=npl+1;
    XB=[];
    YB=[];
    XO=[];
    YO=[];
    [depth,xc,yc,XB,YB,XO,YO]=grid_plane_zeyu(as(ias),Lp,Wp,XB,YB,...
                                               XO,YO,dL,dW,nL,npl);
    k=sum(tSm(1:npl))*NT;
    u=US(k+1:k+sum(tSm(npl+1:npl+1))*NT);
    [U]=split_U(u,Mode)/fact;
    if Mode(3)==0, U(3,:)=U(1,:); end
     xr=[];
     yr=[];
     xmi=Inf;
     xma=-Inf;
     ymi=Inf;
     yma=-Inf;
     k=0;
     for jc=1:Wp(npl)
      for ic=1:nL(npl,jc)
       k=k+1;
       xb=XB((k-1)*5+1:(k-1)*5+4);
       yb=YB((k-1)*5+1:(k-1)*5+4);
%       dt=depth(ic,jc)-dW(npl,jc)*sin(asO(6));
%       zb=[dt depth(ic,jc) depth(ic,jc) dt];
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
       if U(1,k) > 0, U(1,k)=-U(1,k); end;
        if Mode(3) == 0
         fill3(xb,yb,zb,-U(1,k),'EdgeColor','k'); hold on  %
        else
         fill3(xb,yb,zb,U(3,k),'EdgeColor','k'); hold on  %
        end
       plot3(xc,yc,-depth(ic,jc),'ko'), hold on
      end
     end
%     cmax=max([cmax reshape(U1,1,length(U1(:)))]);
%     cmin=min([cmin reshape(U1,1,length(U1(:)))]);
%     u=sqrt(u1.^2+u2.^2+u3.^2); cmax=max(u(:));
     ax1=gca;
     caxis([cmin,cmax]); hold on
%       quiver3([xv(j)],[yv(j)],[zv(j)],[u1(j)/SS],[u2(j)/SS],0,'k'); hold on
     for j=1:length(xv(:))
%      if u(j)>cmax/2
       quiver3([xv(j)],[yv(j)],[zv(j)],[u1(j)/SS],[u2(j)/SS],0,'k'); hold on
%      else
%       quiver3([xv(j)],[yv(j)],[zv(j)],[u1(j)/SS],[u2(j)/SS],0,'w'); hold on
%      end
     end
%     eval(['print ' int2str(i) ' -dpsc2'])

%   text(xc(nL(npl,1),1)-2,yc(nL(npl,1),1), ...
%       -depth(nL(npl,1),1)-3,int2str(i),'color','w');  % plot plane #s
end  % # of sources


axis('equal'), hold on
set(gca,'zlim',[-11 0]);
eval(gotodir0);

return

xb = [-15 -12 -12 -15];
yb = [27 34 34 27];
zb = [0 0 12 12];
fill3(xb,yb,-zb,0,'EdgeColor','k'); hold on  %

axis('equal'), hold on

return

load Hauksson_ridgecrest_17_July.reloc;
x=Hauksson_ridgecrest_17_July;
figure
return

xb = [-15 -12 -12 -15];
yb = [27 34 34 27];
zb = [0 0 12 12];
fill3(xb,yb,-zb,0,'EdgeColor','k'); hold on  %

axis('equal'), hold on

return

load Hauksson_ridgecrest_17_July.reloc;
x=Hauksson_ridgecrest_17_July;
figure
X=x(:,3);
Y=x(:,2);
Z=x(:,4);
%plot(X, Y,'.k'), hold on

d=x(:,13);
h=x(:,14);
m=x(:,15);

  event = find(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84));
%event = find(~(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84)));

%event = find(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84) & d < 5 & (d == 5 & h < 20) );

%plot(X(event), Y(event),'.k'), hold on
%plot(X(event), Y(event),'.k'), hold on
plot3(X(event), Y(event), -Z(event),'.k'), hold on

load Ridgecrest1.trace;
s=Ridgecrest1;

%plot(s(:,1),s(:,2),'or'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

load Ridgecrest2.trace;
s=Ridgecrest2;

%plot(s(:,1),s(:,2),'ob'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

saveas(gcf,'seism.fig')
return

x1(:,1) =X(event);  
x1(:,2) =Y(event);

%save eq2.dat x1 -ascii
