% Plots output for the HM model
% LOAD INPUTS
clear all;
[dpath]=read_line('/Users/fialko/w/matlab/work/pathname.dat');
%dpath='/utah/usr/fialko/insar/hector/model/layered/';
gotodir=['cd ' dpath];
eval(gotodir);
rlks=8;
fid = fopen([dpath 'looks.in'],'wt'); % define # of looks
fprintf(fid,'%d \n',rlks);
status = fclose(fid);
shift=89; %define the InSAR ramp shape
fid = fopen([dpath 'shift.in'],'wt'); % and write it to a file
fprintf(fid,'%d \n',shift);
status = fclose(fid);
BOX=0;            % 1 - use data from a limited area (specified below)
fid = fopen([dpath 'data_box.in'],'wt'); % and write it to a file
%if BOX==1, fprintf(fid,'%f %f %f %f \n',-15,10,-20,15); end;
%if BOX==1, fprintf(fid,'%f %f %f %f \n',-116.6 -116.3 34.4 34.8); end;
if BOX==1, fprintf(fid,'%f %f %f %f \n',-117, -116, 39.8, 35.2); end;
status = fclose(fid);

load_data_2plot;
%return
load([dpath 'orig.ll']);
xo=orig(1); yo=orig(2); 
[X,Y]=utm2ll(xo,yo,0,1);
xo=X; yo=Y;
[X1,Y1]=utm2ll(xG,yG,0,1);
xG=(X1-X)*1e-3;
yG=(Y1-Y)*1e-3;

[X1,Y1]=utm2ll(cont(:,1),cont(:,2),0,1);
cont(:,1)=(X1-X)*1e-3;
cont(:,2)=(Y1-Y)*1e-3;

[X1,Y1]=utm2ll(camp(:,1),camp(:,2),0,1);
camp(:,1)=(X1-X)*1e-3;
camp(:,2)=(Y1-Y)*1e-3;

dat_ph=0;
dat_az=0;

fault_tr=1; % put fault_tr=1 if want to plot the geologic fault trace
pr=1;       % put pr=1 if want to print output into files
wr=0;       % put wr=1 if want to write residual 
units=0;
%[xo,yo]=utm2ll(-116.27,34.595,0,1); %epicenter coordinates

%load MODEL3_L.mat;
% Load model parameters
clear US as Lp Wp dL nL dW tSm Mode MI shift;

load MODEL3_Z.mat;
MM=8;
shift=89;
fid = fopen([dpath 'shift.in'],'wt'); % and write it to a file
fprintf(fid,'%d \n',shift);
status = fclose(fid);

%global MA NS NM Mode src_type MI tSm shift FA TP GF edcmp
eps=1e-10;              % small number


%x1=-40:5:40; y1=-20:5:60;
%[xx,yy]=meshgrid(x1,y1);
%xG=xx(:); yG=yy(:);

[G1,G2]=size(xG);
uztot=zeros(size(xG));
range_sum=zeros(size(xd));
ranged=[];
%plot initial model
 XB=[];
 YB=[];
 XO=[];
 YO=[];
npl=0;
MM=0;

if fault_tr==1
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
% [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');
 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
% [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');

end


bgrey = [1 1 1]*0.7;
cbar = 'v';

Fontsize=11;

%return
NT=sum(Mode);
matrl=[1 1 0.25];

fact = [0.5 0.6 0.7 0.8 0.9  1 ];



%for it=1:length(fact)
for it=6:6
fig1 = figure(1);
clf

cmax=-Inf;
cmin=Inf;
% Do horizontals:
MM=0;
uxtot=zeros(size(xG));
uytot=zeros(size(xG));
uztot=zeros(size(xG));
tpG=zeros(size(xG));

Mo=0;
Moss=0;

%plot initial model
XB=[];
YB=[];
XO=[];
YO=[];
npl=0;
k=0;
strace=[];
num_seg=0;
for i= 1:NS
%for i= 1:1
 if i > 1, MM=MI(i-1); end
  ias=(1:MI(i))+MM*(i-1);
  if src_type(i) == 'Y' | src_type(i) == 'y' 
    fprintf('huilo\n');
    [ranged,ux,uy,uz]=fcn_yang(as(ias),x,y,matrl,s_uv,tp);
    uxtot=uxtot+ux;
    uytot=uytot+uy;
  elseif src_type(i) == 'O' | src_type(i) == 'o' 
    asO=as(ias);
    strike=asO(7)/pi*180;
    asO=as(ias);
    npl=npl+1;
    [depth,xc,yc,XB,YB,XO,YO]=grid_plane_zeyu1(asO,Lp,Wp,XB,YB,XO,YO,dL,dW,nL,npl);

%depth(1,1)
    if depth(1,1)<0.5
     x1=xc(1,1)+asO(4)/2*cos(pi/2-asO(7));
     x2=xc(1,1)-asO(4)/2*cos(pi/2-asO(7));
     y1=yc(1,1)+asO(4)/2*sin(pi/2-asO(7));
     y2=yc(1,1)-asO(4)/2*sin(pi/2-asO(7));
     strace=[strace x1 x2 y1 y2];
     num_seg=num_seg+1;
    end

    if depth(1,1)<0.5
%     plot(xc(1,1),yc(1,1),'ok'), hold on
%      line(xc(1,1),yc(1,1),'LineStyle','none','Marker','o','MarkerSize',10,'MarkerFaceColor','k'), hold on
    end
    k=sum(tSm(1:npl))*NT;
    u=US(k+1:k+sum(tSm(npl+1:npl+1))*NT);
    [U]=split_U(u,Mode);
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
       xr=xG-xc(ic,jc);
       yr=yG-yc(ic,jc);
       for j=1:3

slip=0;
       if Mode(j)==1
        k=k+1;
%if strike > 300 & j==1, slip=-4e1; end  % equiv. M6.05
%if strike > 300 & j==1, slip=-2e1; end  % equiv. M6.05
if strike > 300, slip=0e1; end  % equiv. M6.0
%if strike < 300 & j==1, slip=1e1; end
%if strike > 300, slip=U(j,1)*1e-1; end
if strike < 300, slip=U(j,1)*fact(it); end
        [ux,uy,uz]=calc_green(asO(8),slip,xr,yr,matrl(3),asO(6),...
                              depth(ic,jc),dL(npl,jc),dW(npl,jc),j,asO(7),tpG);
        if j==1
         Uss=slip;
         slip=0;
        end
        if j==2
         Uds=slip;
         slip=0;
        end
        uxtot=uxtot+ux;
        uytot=uytot+uy;
        uztot=uztot+uz;
       end
      end
      Moss=Moss+Uss*dL(npl,jc)*dW(npl,jc);
      Mo=Mo+sqrt(Uss^2+Uds^2)*dL(npl,jc)*dW(npl,jc);
     end
    end
  end
end

Mw(it)=(log10(Mo*3.3e14)-9.1)/1.5;
%fprintf('\n strike-slip seismic moment= %e \n',Moss*3.3e14);
%fprintf('\n Total seismic moment= %e \n',Mo*3.3e14);
%fprintf('\n Total moment magnitude= %e \n',Mw(it));

run_time=cputime;
%pause


% npl=0;
% strace=[];
% for i= 1:NS
%  if i > 1, MM=MI(i-1); end
%  ias=(1:MI(i))+MM*(i-1);
%  if src_type(i) == 'O' | src_type(i) == 'o' 
%    asO=as(ias);
%    x1=asO(1)+asO(4)/2*cos(pi/2-asO(7));
%    x2=asO(1)-asO(4)/2*cos(pi/2-asO(7));
%    y1=asO(2)+asO(4)/2*sin(pi/2-asO(7));
%    y2=asO(2)-asO(4)/2*sin(pi/2-asO(7));
%    if units==1
%     [X1,Y1]=utm2ll(x1*1e3+xo,y1*1e3+yo,UTM,2);
%     x1=X1;
%     y1=Y1;
%     [X1,Y1]=utm2ll(x2*1e3+xo,y2*1e3+yo,UTM,2);
%     x2=X1;
%     y2=Y1;
%    end  
%    strace=[strace x1 x2 y1 y2];
%    npl=npl+1;
%  end
% end

xlim=[-50 45]; ylim=[-30 65];
clim=[-5 5];
save SOL_GPS xG yG uxtot uytot uztot

%return

%set(fig1,'Name','Camp. GPS data & predictions',...
%    'PaperPosition',[0.25 1.5 8 8],...
%    'Position',[10 20 660 600],'number','off');
if fault_tr==1
 [i]=plot_faults(faults,dim,xlim,ylim,'m');
end

for j=1:num_seg
  xb=strace((j-1)*4+1:(j-1)*4+2);
  yb=strace((j-1)*4+3:(j-1)*4+4);
  line(xb,yb,'Color',bgrey,'LineStyle','-','LineWidth',3),hold on
end

%return

%A=0.8;
%A=0.4;
%hscale=[20]*A;
%htext = ['20 cm'];
%b=-30;
%xe=-50; ye=90;
%xlim=[-80 70]; ylim=[-90 60];
%A=0.2;  % co-seis

A=3;
hscale=[5]*A;
htext = ['50 mm'];
b=0;
xe=33; ye=-22;
%q1=quiver([-10]+xe,[55]+ye,hscale,[0],0,'k'); hold on
text(-7+xe,58+ye,htext,'Color','k','FontSize',14), hold on
V(1)=hscale;
V(2)=[0];
quiver([-10]+xe,[55]+ye,V(1),V(2),0,'k','AutoScale','off','MaxHeadSize',5/norm(V)); hold on


%quiver(cont(:,1),cont(:,2),A*cont(:,3),A*cont(:,4),0, 'r'), hold on
%quiver(camp(:,1),camp(:,2),A*camp(:,3),A*camp(:,4),0, 'b'), hold on

insd=find(xG>xlim(1) & xG<xlim(2) & yG>ylim(1) & yG<ylim(2)); 
for i=1:length(insd(:))
  x5=[xG(insd(i)) xG(insd(i))+A*uxtot(insd(i))];
  y5=[yG(insd(i)) yG(insd(i))+A*uytot(insd(i))];
% line(x5,y5,'LineStyle','-','Color','k','Marker','.'), hold on
% plot(xG(i)+A*uxtot(i),yG(i)+A*uytot(i),'kp'), hold on
%quiver_arrow(xG,yG,A*uxtot,A*uytot,'k'), hold on
end
%quiver_arrow(xG(insd),yG(insd),A*uxtot(insd),A*uytot(insd),'k'), hold on

% quiver(xG,yG,A*cont(:,3)*1e2,A*cont(:,4)*1e2,'r'), hold on
numpo=91;
for i=1:m1
 xf=cont(i,1)+A*cont(i,3);
 yf=cont(i,2)+A*cont(i,4);
 x5=[cont(i,1) xf];
 y5=[cont(i,2) yf];
 Dxel=3*A*cont(i,5)/(numpo-1);
 xel=[(xf-A*cont(i,5)):Dxel:(xf+A*cont(i,5))];
 yel=real(A*cont(i,6)*sqrt(1-((xel-xf)/A/cont(i,5)).^2));
 yel=[-yel yel]+yf;
 xel=[xel xel(length(xel):-1:1)];
 line(x5,y5,'LineStyle','-','Color','b','Marker','.'), hold on
 plot(cont(i,1)+A*cont(i,3),cont(i,2)+A*cont(i,4),'b*'), hold on
 line(xel,yel,'LineStyle','-','Color','b'), hold on
end
%return

for i=1:m3
 xf=camp(i,1)+A*camp(i,3);
 yf=camp(i,2)+A*camp(i,4);
 x5=[camp(i,1) xf];
 y5=[camp(i,2) yf];
 Dxel=3*A*camp(i,5)/(numpo-1);
 xel=[(xf-A*camp(i,5)):Dxel:(xf+A*camp(i,5))];
 yel=real(A*camp(i,6)*sqrt(1-((xel-xf)/A/camp(i,5)).^2));
 yel=[-yel yel]+yf;
 xel=[xel xel(length(xel):-1:1)];
 x5=[camp(i,1) camp(i,1)+A*camp(i,3)];
 y5=[camp(i,2) camp(i,2)+A*camp(i,4)];
 line(x5,y5,'LineStyle','-','Color','k','Marker','.'), hold on
 plot(camp(i,1)+A*camp(i,3),camp(i,2)+A*camp(i,4),'kp'), hold on
 line(xel,yel,'LineStyle','-','Color','k'), hold on
end

xe = -23; ye = 95;
x5=[10 15]+xe;

a=-36+ye;y5=[a a];
line(x5,y5,'LineStyle','-','Color','b','Marker','.'), hold on
plot(x5(2),y5(2),'b*'), hold on
text(x5(2)+5,y5(2),'continuous GNSS','Color','b','Fontsize',12);

a=-39+ye;y5=[a a];
line(x5,y5,'LineStyle','-','Color','k','Marker','.'), hold on
plot(x5(2),y5(2),'kp'), hold on
text(x5(2)+5,y5(2),'campaign GNSS', ...
     'Color','k','Fontsize',12);


a=-42+ye;y5=[a a];
%%line(x5,y5,'LineStyle','-','Color','k','Marker','.'), hold on
%%plot(x5(2),y5(2),'kp'), hold on
%quiver_arrow([x5(1)],[y5(1)]-3,[x5(2)-x5(1)],[y5(2)-y5(1)],'r'), hold on
V(1)=[x5(2)-x5(1)]+1;
V(2)=[y5(2)-y5(1)];
q1=quiver([x5(1)],[y5(1)],V(1),V(2),0,'r','AutoScale','off','MaxHeadSize',10/norm(V)); hold on
text(x5(2)+5,y5(2),'Model','Color','r','Fontsize',12);
%q.MaxHeadSize = 0.1;

q=quiver(xG(insd),yG(insd),A*uxtot(insd),A*uytot(insd),0,'r'); hold on
q.MaxHeadSize = 0.05;

text(9,14,'PNCL','Color','k','Fontsize',12);

%plot(0,0,'rp','MarkerSize',10,'MarkerFaceColor','r'), hold on  % epicenter location
axis('equal')
xlabel('Eastings, km','Fontsize',14);
ylabel('Northings, km','Fontsize',14);
set(gca,'Fontsize',14,'XLim',xlim,'YLim',ylim,'box','on');

orientation = [0.05, 0.1, 7, 7; 0.1, 0.1,7, 7];
computerNames = {'foxbat'; 'flagon'};
[~, thiscomputer] = system('hostname');
option = strcmp(deblank(thiscomputer), computerNames);
set(gcf,'defaulttextinterpreter','none','Visible','off','PaperOrientation', 'landscape','PaperPositionMode','manual', 'PaperPosition',orientation(option, :))


if pr==1
% print -r300 % set resolution
 print('-dpng','-r300','GPS_M6')
 % print GPS -dpsc2
 %print GPS -dill
% saveas(gcf,'GPS_M6.png');
% saveas(gcf,'junk.png');
end

allgps=[cont' camp']';

misfit(it) = sum((allgps(:,3)-uxtot).^2 +(allgps(:,4)-uytot).^2);
%[fact(it) misfit(it) Mw(it)]

%fprintf('%e %e \n',fact(it),Mw(it));
fprintf('%e %e %e \n',fact(it),misfit(it),Mw(it));


end

