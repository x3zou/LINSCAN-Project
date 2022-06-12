clear all;
close all;

bgrey = [1 1 1];
wid=2;
fault_tr = 1; % flag for plotting the fault trace

write_line('/Users/fialko/matlab/work/pathname.dat','/Users/fialko/w/tex/papers/Ridgecrest1/figs1/');
[dpath]=read_line('/Users/fialko/matlab/work/pathname.dat');
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

M_T=zeros(1,6);
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
U(1,k)=1;
U(2,k)=0;
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
%fill3(xb,yb,zb,U(i),'EdgeColor','k','FaceAlpha',0.4); hold on  %
%         fill3(xb,yb,zb,-U(1,k),'EdgeColor','k'); hold on  %
%        else
%         fill3(xb,yb,zb,U(3,k),'EdgeColor','k'); hold on  %
%        end
       
ut=sqrt(U(1,k)^2+U(2,k)^2);

strike = asO(7)/pi*180;

%if strike > 280
%if strike < 80
mo=dW(npl,jc)*dL(npl,jc)*ut;
rake=atan2(U(2,k),U(1,k)); % for now
[momt]=mt(asO(7),asO(6),rake,mo);
M_T=M_T+momt;

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

fprintf('\n %f %f %f %f %f %f \n',M_T(1:6));
axis('equal'), hold on
set(gca,'zlim',[-11 0]);
eval(gotodir0);

M_T

