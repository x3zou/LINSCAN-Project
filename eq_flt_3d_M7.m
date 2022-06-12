clear all;
close all;

load Hauksson_ridgecrest_17_July.reloc;
x=Hauksson_ridgecrest_17_July;
X=x(:,3);
Y=x(:,2);
Z=x(:,4);
%plot(X, Y,'.k'), hold on
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

event = find(X> -117.9 & X < -117.2 & Y < 35.85 & Y > 35.75);
x=X(event); 
y=Y(event);
[X1,Y1]=utm2ll(x,y,0,1);
xU=(X1-xo)*1e-3;
yU=(Y1-yo)*1e-3;
 
%plot(xU, yU,'.k'), hold on
%return

plot3(xU, yU, -Z(event),'.k'), hold on

load Ridgecrest1.trace;
x=Ridgecrest1;
X=x(:,1);
Y=x(:,2);
Z=0*x(:,1);
event = find(X> -117.9 & X < -117.2 & Y < 35.85 & Y > 35.75);

x=X(event);
y=Y(event);

[X1,Y1]=utm2ll(x,y,0,1);
xU=(X1-xo)*1e-3;
yU=(Y1-yo)*1e-3;

plot3(xU, yU, Z(event),'or'), hold on

axis('equal')
%return

%write_line('/Users/fialko/matlab/work/pathname.dat','/Users/fialko/w/tex/papers/Ridgecrest1/model/');
write_line('/Users/fialko/w/matlab/work/pathname.dat','/Users/fialko/w/tex/papers/Ridgecrest1/figs1/');
[dpath]=read_line('/Users/fialko/w/matlab/work/pathname.dat');
gotodir=['cd ' dpath];
eval(gotodir);
gotodir0=['cd /Users/fialko/w/tex/papers/Ridgecrest1/figs1'];

global MA NS NM Mode src_type MI tSm shift FA TP GF edcmp
load_param_inv; % load the initial model
topmin=0;
MW=0; AW=0; CW=0; GW=0;
setup_inp_inv;  % prepare inputs for the inversion
%eps=1e-10;              % small number
%pack;
%return

save MODEL as MI NS shift Mode HL GF TP;
%return

load ../figs/MODEL3_Z.mat;

return
options = optimset('Display','final');
ndat=1; mdat=1;  % these are some dummy variables
    fpar(1)=dl;
    fpar(2)=dw;

eval(gotodir);

dW=[]; dL=[]; nL=[]; Wp=[]; Lp=[]; 
MM=0; npl=0; 
% divide faults into subfaults
%NS=6;
for i= 1:NS 
 if i > 1, MM=MI(i-1); end
 ias=(1:MI(i))+MM*(i-1);
 if src_type(i) == 'O' | src_type(i) == 'o' 
   asO=as(ias);
   npl=npl+1;
%   if i==7, dw=dW(1)*FA; dl=dw; end 
   [nL,dW,dL,Wp,Lp,tSm(i+1)]=expl_plane(dL,nL,dW,asO,Wp,Lp,dw,dl,npl);
 end
end
NT=sum(Mode);

npl=0;                 % plot slip distribution
SS=2;                 % slip scale 
ut1=[];
ut2=[];
MM=0;
xs=[];xv=[];yv=[];zv=[];u1=[];u2=[];u3=[];

U(1)=0; U(2)=1;
for i= 1:NS
 if i > 1, MM=MI(i-1); end
  ias=(1:MI(i))+MM*(i-1);
    asO=as(ias);
    npl=npl+1;
    XB=[];
    YB=[];
    XO=[];
    YO=[];
    [depth,xc,yc,XB,YB,XO,YO]=grid_plane(as(ias),Lp,Wp,XB,YB,...
					       XO,YO,dL,dW,nL,npl);
    k=sum(tSm(1:npl))*NT;
%    u=US(k+1:k+sum(tSm(npl+1:npl+1))*NT);
%    [U]=split_U(u,Mode)/fact;
%    if Mode(3)==0, U(3,:)=U(1,:); end  
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
       dt=-depth(ic,jc)-dW(npl,jc)*sin(asO(6));
       zb=[dt -depth(ic,jc) -depth(ic,jc) dt];
       xv=[xv 0.5*(xb(1)+xb(3))];
       yv=[yv 0.5*(yb(1)+yb(3))];
       zv=[zv 0.5*(dt-depth(ic,jc))];
%       x1=U(2,k)*cos(asO(6)); y1=U(1,k);
%       x2=x1*cos(asO(7))+y1*sin(asO(7));
%       y2=-x1*sin(asO(7))+y1*cos(asO(7));
%       u1=[u1 x2];
%       u2=[u2 y2];
%       u3=[u3 U(2,k)*sin(asO(6))];
%        if Mode(3) == 0 
fill3(xb,yb,zb,U(i),'EdgeColor','k','FaceAlpha',0.4); hold on  %
%         fill3(xb,yb,zb,-U(1,k),'EdgeColor','k'); hold on  %
%        else
%         fill3(xb,yb,zb,U(3,k),'EdgeColor','k'); hold on  %
%        end
      end
     end
%     cmax=max([cmax reshape(U1,1,length(U1(:)))]);
%     cmin=min([cmin reshape(U1,1,length(U1(:)))]);
%     u=sqrt(u1.^2+u2.^2+u3.^2); cmax=max(u(:));
     ax1=gca;
%     caxis([cmin,cmax]); hold on
%     eval(['print ' int2str(i) ' -dpsc2'])

%   text(xc(nL(npl,1),1)-2,yc(nL(npl,1),1), ...
%	-depth(nL(npl,1),1)-3,int2str(i),'color','w');  % plot plane #s
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
