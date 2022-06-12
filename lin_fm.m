clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
%minlon = -32; maxlon = -25; minlat = 32; maxlat = 40; 
% %minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
% %minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
Fs=14;

deg2rad=pi/180.;

for i=1:0
%s=340;d=125;r=85;
s=rand(1)*360; d=-90+rand(1)*180; r=-180+rand(1)*360; 
% s=87; d=46; r=45;
% s=84; d=35; r=54;
%s=309; d=60; r=145;
%s=190.925; d=42.4899; r=-20.9735; 
figure(2)
clf
subplot(1,3,1)
bb1([s d r], 0, 0, 1.5, 0, 'b'), hold on
axis equal;
set(gca,'xlim',[-1.7 1.7],'ylim',[-1.7 1.7],'box','on','FontSize',Fs);
text(-1,3,num2str(s)), hold on
text(-1,2.8,num2str(d)), hold on
text(-1,2.6,num2str(r)), hold on

[fm]=sdr2mij([s d r],deg2rad);
%figure(2)
subplot(1,3,2)

bb1(fm, 0, 0, 1.5, 0, 'b'), hold on
axis equal;
set(gca,'xlim',[-1.7 1.7],'ylim',[-1.7 1.7],'box','on','FontSize',Fs);
 
subplot(1,3,3)
focalmech(fm, 0, 0, 3.2, 'r')
axis equal;
set(gca,'xlim',[-1.7 1.7],'ylim',[-1.7 1.7],'box','on','FontSize',Fs);

% %function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
% [s1,d1,r1,s2,d2,r2] = mij2sdr2(fm);
% %size(fm)
% %fm
% %[s1,d1,r1] = mij2sdr1(fm);
% %figure(3)
% subplot(1,3,3)
% bb([s1 d1 r1], 0, 0, 1.5, 0, 'b'), hold on
% text(-1,3,num2str(s1)), hold on
% text(-1,2.8,num2str(d1)), hold on
% text(-1,2.6,num2str(r1)), hold on
% text(1,3,num2str(s2)), hold on
% text(1,2.8,num2str(d2)), hold on
% text(1,2.6,num2str(r2)), hold on
%fprintf('%f,%f,%f,%f,%f,%f\n',fm(:));
%fprintf('%f,%f,%f,%f,%f,%f\n',fm(),fm(),fm(),fm(),fm(),fm());
%fprintf('%f,%f,%f,%f,%f,%f\n',fm(2),fm(3),fm(1),-fm(6),fm(4),-fm(5));

%bb([200 90 0], 0, 0, 1.5, 0, 'b'), hold on
%   axis equal;
pause
end
%return

pr=1;
fault_tr = 0; % flag for plotting the fault trace
% regional
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
grey = [1 1 1]*0.6;
% %grey = [0.4,0.4,0.4];
% 
%load lin.utm;
%X=lin;

f = 'sc2019_hash_ABCD_so.focmec.scedc.txt';
% f = 'sc.focmec.scedc.txt';
fid = fopen(f);
C = textscan(fid,'%d %d %d %d %d %f %d %f %f %f %f %f %f %f %f %f %f %f %d %f %s', 'delimiter', '\n');
fclose(fid);

wid=length(C(1,:));
N=cell2mat(C(:,8:10));    % coordinates
FM=cell2mat(C(:,11:14)); % focal mech (M strike dip rake)
yr=cell2mat(C(:,1));     % year
mo=cell2mat(C(:,2));     % month
L=C{:,wid};              % quality flag

%ind=find((yr< 2019 |(yr==2019 & mo<7)) );
ind=find(yr< 2020);
%ind=find(strtrim(char(L)) == 'A');
%ind=find(strtrim(char(L)) == 'A' | strtrim(char(L)) == 'B');  % select high-quality data

x=N(ind,2);
y=N(ind,1);
z=-N(ind,3);
FM=FM(ind,1:4);
qn=0*x;
Q=strtrim(char(L(ind)));     
% assign numerical score to quality flag 
ind=find(Q == 'A');
qn(ind)=4;
ind=find(Q == 'B');
qn(ind)=3;
ind=find(Q == 'C');
qn(ind)=2;
ind=find(Q == 'D');
qn(ind)=1;

%return


[xx,yy]=utm2ll(x(:),y(:),0,1);
% !!!!!
xshift=0.0; % km; stress drop locations are systematically to the E from Lin relocated catalog

xu=reshape((xx-xo)*1e-3,size(x)) - xshift;
yu=reshape((yy-yo)*1e-3,size(x));

ind=find(xu < maxlon & xu > minlon & yu < maxlat & yu > minlat );

length(ind)
%#line(xu(ind),yu(ind),'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor',grey,'MarkerEdgeColor',grey), hold on

xlin=[xu yu]; 
matr = [xu(ind) yu(ind)  z(ind) ];
save jj matr -ascii;
%return

D = 0.5; % distance threshold for selecting nearby points
style = ':';
color = 'k';
width = 1;
color = [rand(1) rand(1) rand(1)];
X=[]; Y=[]; Len=[]; S=[];

if fault_tr==1

 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');

end

f = 'lin_segments.dat';
%f = 'lin_segments.dat2';
if isfile(f)
 fid = fopen(f,'r');
 load(f);
 [m,n]=size(lin_segments);
 x=lin_segments;
 X=0.5*(x(:,1)+x(:,3));
 Y=0.5*(x(:,2)+x(:,4));
 Len=sqrt((x(:,4)-x(:,2)).^2+(x(:,3)-x(:,1)).^2);
 s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
 ind=find(s>0 & s<90);
 s(ind)=90-s(ind);
 ind=find(s>90);
 s(ind)=450-s(ind);
 ind=find(s<0 & s>-90);
 s(ind)=270-s(ind);
 ind=find(s<-90);
 s(ind)=-90-s(ind);
 S=s;

%  for i=1:m
%   line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
%  end
 fclose(fid);
end

%pause

fid=fopen('fc.dat','wt');

% % plot all of the clusters:
x0=X; y0=Y; L=Len; strike=S;
%for i=1:0
for i=1:length(x0)
%for i=length(x0):-1:1
 color = [rand(1) rand(1) rand(1)];  
 color = [0 0 0];  
 color = 'k';  
 width=1;
 %figure(1)
 %clf
 plot_fm(i,x0(i),y0(i),xlin,FM,strike(i),L(i),D,color,style,width,fid,qn);
 %tit=[cellstr(['segment number : ' num2str(i)])];
 %title(char(tit),'Fontsize',14);
 %axis equal;
 %pause
end

fclose(fid);


axis('equal');
xlabel('Eastings, km');
ylabel('Nortings, km');
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);

%saveas(gcf, 'lin_fm', 'png')
print(gcf,'-dpng','-r300','lin_fm');
return

% text(-45,5,'1984-2011','FontSize',Fs);
load preeq_ross.utm;
x=preeq_ross;
%line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',0.1,'MarkerFaceColor','m','MarkerEdgeColor','b'), hold on

axis('equal');
xlabel('Eastings, km');
ylabel('Nortings, km');
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
%legend('earthquakes','\sigma_1','$\dot{\epsilon}$', 'Interpreter','latex');

saveas(gcf, 'lin_trends', 'png')
return


function [x1,y1,x2,y2] = plot_line(x,y,strike,L,color,style,width,fid)
 % Compute end points
  x1 = x + L/2*sind(strike);
  x2 = x - L/2*sind(strike);
  y1 = y + L/2*cosd(strike);
  y2 = y - L/2*cosd(strike);
  if fid > 0
  fprintf(fid,'%f %f \n',x1,y1);
  fprintf(fid,'%f %f \n',x2,y2);
  fprintf(fid,'> \n',color);
  end
%  fprintf(fid, '%s \n',color);


 % line([x1 x2],[y1 y2],'lineWidth',width,'lineStyle',style,'Color',color), hold on
%  line([x1 x2],[y1 y2],'lineWidth',1,'lineStyle','--','Color','m'), hold on
end



function [ind,d] = point_to_line(pt, v1, v2, D)
% function [ind,d] = point_to_line(pt, v1, v2, D)
% return points within distance D to a finite line segment

 L = norm(v1 - v2);

 v1 = repmat(v1,size(pt,1),1); % 1st end of line point
 v2 = repmat(v2,size(pt,1),1); % 2nd end of line point
 a = v1 - v2;
 b = pt - v2;
 c = pt - v1;
 d1= sqrt(sum(c.*c,2)); % square of distance to one end of the line segment
 d2= sqrt(sum(b.*b,2)); % square of distance to the other end of the line segment

% d = sqrt(sum(cross(a,b,2).^2,2) ./ sum(a.^2,2)); %distance to line
 d = sqrt(sum(cross(a,b,2).^2,2))/L; %distance to line

 ind=find(d < D & d1+d2 < L+2*D); % index to nearby points

end


function [str,dip,rake] = mij2sdr1(mxx,myy,mzz,mxy,mxz,myz)
%function [str,dip,rake] = mij2sdr1(mxx,myy,mzz,mxy,mxz,myz)
%
%   INPUT
%	mij - siz independent components of the moment tensor
%
%   OUTPUT
%	str - strike of first focal plane (degrees)
%	dip - dip of first focal plane (degrees)
%	rake - rake of first focal plane (degrees)
%
%Adapted from code, mij2d.f, created by Chen Ji and given to me by Gaven Hayes.
%

a = [mxx mxy mxz; mxy myy myz; mxz myz mzz];
[V,d] = eig(a);

D = [d(3,3) d(1,1) d(2,2)];
V(2:3,1:3) = -V(2:3,1:3);
V = [V(2,3) V(2,1) V(2,2); V(3,3) V(3,1) V(3,2); V(1,3) V(1,1) V(1,2)];

IMAX = find(D == max(D));
IMIN = find(D == min(D));
AE = (V(:,IMAX)+V(:,IMIN))/sqrt(2.0);
AN = (V(:,IMAX)-V(:,IMIN))/sqrt(2.0);
AER = sqrt(AE(1)^2+AE(2)^2+AE(3)^2);
ANR = sqrt(AN(1)^2+AN(2)^2+AN(3)^2);
AE = AE/AER;
AN = AN/ANR;
if (AN(3) <= 0.)
	AN1 = AN;
	AE1 = AE;
else
	AN1 = -AN;
	AE1 = -AE;
end 
[ft,fd,fl] = TDL(AN1,AE1);
str = 360 - ft;
dip = fd;
rake = 180 - fl;

end 


function [FT,FD,FL] = TDL(AN,BN)
XN=AN(1);
YN=AN(2);
ZN=AN(3);
XE=BN(1);
YE=BN(2);
ZE=BN(3);
AAA=1.0E-06;
CON=57.2957795;
if (abs(ZN) < AAA)
	FD=90.;
	AXN=abs(XN);
	if (AXN > 1.0) 
		AXN=1.0;
	end
	FT=asin(AXN)*CON;
	ST=-XN;
	CT=YN;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	FL=asin(abs(ZE))*CON;
	SL=-ZE;
	if (abs(XN) < AAA) 
		CL=XE/YN;
	else
		CL=-YE/XN;
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
else
	if (-ZN > 1.0) 
		ZN=-1.0;
	end
	FDH=acos(-ZN);
	FD=FDH*CON;
	SD=sin(FDH);
	if  (SD == 0)
		return;
	end 
	ST=-XN/SD;
	CT=YN/SD;
	SX=abs(ST);
	if (SX > 1.0) 
		SX=1.0;
	end
	FT=asin(SX)*CON;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	SL=-ZE/SD;
	SX=abs(SL);
	if (SX > 1.0) 
		SX=1.0;
	end
	FL=asin(SX)*CON;
	if (ST == 0) 
		CL=XE/CT;
	else
		XXX=YN*ZN*ZE/SD/SD+YE;
		CL=-SD*XXX/XN;
		if (CT == 0) 
			CL=YE/ST;
		end
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
end 

end

function plot_fm(k,x0,y0,xlin,FM,strike,L,D,color,style,width,fid,qn)

%[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,fid);
[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,'k',style,width,0);
%plot([x1 x2],[y1 y2],'ko'), hold on
v1=[x1 y1 0];  v2=[x2 y2 0]; 

nx=[]; ny=[];
pts=[xlin(:,1) xlin(:,2) 0*xlin(:,1)];

[ind] = point_to_line(pts, v1, v2, D);
nx = xlin(ind,1);
ny = xlin(ind,2);
Q = mean(qn(ind));

%#line(nx,ny,'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor',color,'MarkerEdgeColor',color), hold on

%[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,'k',style,4,fid);

fm=FM(ind,1:4);

arr=[nx ny fm];

deg2rad=pi/180.;
comp=zeros(1,6);
sm=[];
num_events=length(arr(:,1));
 for i=1:num_events
  [fm]=sdr2mij(arr(i,4:6),deg2rad);
  Mo=10^(1.5*(arr(i,3)+10.7));
  fm=fm*Mo;
%  Mo=sqrt(sum(fm(:).^2)/2);
%  fm=fm/Mo;
  comp=comp+fm;
  % this is to estimate variability:
  [s1,d1,s2,d2] = mij2sd(fm);
  sm=[sm min([s1 s2])];
 end
stdm=std(sm); 

%Mo=max(abs(comp(1:3)));
Mo=sqrt(sum(comp(:).^2)/2);
comp=comp/Mo;
%fprintf('%f,%f,%f,%f,%f,%f\n',comp(:));
%fprintf('%f,%f,%f,%f,%f,%f\n',comp(2),comp(3),comp(1),-comp(6),comp(4),-comp(5));
%comp=comp-sum(comp(1:3));
% Get the best-fit double couple
%[V,D] = eig(comp);
%e=trace(D)/3;
%function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)

M = [comp(1) comp(4) comp(5); comp(4) comp(2) comp(6); comp(5) comp(6) comp(3)];

%return

%sum(comp(1:3))
%pause
 
[s1,d1,r1,s2,d2,r2] = mij2sdr2(comp);

% % average mechanism:
style = '--';
%color = 'b';
width = 3;
%x0=mean(arr(:,1));
%y0=mean(arr(:,2));
%ds=rem(abs(strike-s1),90);
ds=abs(strike-s1);
if abs(rem(round(ds/90),2))==0
%if ds<45
  s=s1;
  r=r1;
else
  s=s2;
  r=r2;
end

% strike from focal mechanism:
%[x1,y1,x2,y2]=plot_line(x0,y0,s,L,color,style,width,fid);
%[x1,y1,x2,y2]=plot_line(x0,y0,s,L,'m',style,width,0);

xm=mean(nx)+3; ym=mean(ny)-3;
%#line([x0 xm],[y0 ym],'Color','k','LineStyle','-','LineWidth',0.5);
%figure(5)
%clf
%bb(comp, mean(nx)+3, mean(ny)-3, 1.5, 0, 'b'), hold on

focalmech(comp, xm, ym, 3.5, 'b','dc'), hold on
%focalmech(comp, xm, ym, 3.5, 'b'), hold on
% title('full');
% figure(6)
% clf
% bb([s1 d1 r1], mean(nx)+3, mean(ny)-3, 1.5, 0, 'b'), hold on
% title('DC');
% pause
% return

%text(mean(nx)+1.5, mean(ny)-1.5,num2str(num_events),'FontSize',12), hold on
% text(mean(nx)+0.5, mean(ny)+0.5,num2str(k),'FontSize',11), hold on
% text(mean(nx)+0.57, mean(ny)+0.57,num2str(Q),'FontSize',11,'Color',color), hold on
% text(mean(nx)+1.2, mean(ny)+0.59,num2str(stdm),'FontSize',14,'Color','k'), hold on
% text(mean(nx)+1.2, mean(ny)+0.39,num2str(mean(sm)),'FontSize',14,'Color','k'), hold on
% text(mean(nx)+1.2, mean(ny)+0.19,num2str(s1),'FontSize',14,'Color','k'), hold on
%text(mean(nx)+0.5, mean(ny)-0.5,num2str(strike),'FontSize',11,'Color',color), hold on
%text(mean(nx)+0.5, mean(ny)-0.53,num2str(s1),'FontSize',11,'Color',color), hold on
%text(mean(nx)+0.5, mean(ny)-0.56,num2str(ds),'FontSize',11,'Color',color), hold on

%[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,0);
[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,'g',style,width,0);

%fprintf('%d %f \n',k,strike);

  if fid > 0
 %  fprintf(fid,'%f %f %f %f \n',x1,y1,x2,y2);
%   fprintf(fid,'%f %f %f %f %f %f %f \n',x1,y1,x2,y2,s1,d1,r1);
%   fprintf(fid,'%f %f %f %f %f %f %f %f %f %f \n',x1,y1,x2,y2,comp(:));
   fprintf(fid,'%8.4f %8.4f %8.4f %8.4f %f %f %f %f %f %f \n',x1,y1,x2,y2,comp(:));
  end

end

