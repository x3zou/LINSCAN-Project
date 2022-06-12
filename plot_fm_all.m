clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
% %minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
% %minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
Fs=14;

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
FM=cell2mat(C(:,12:14)); % focal mech (strike dip rake)
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
FM=FM(ind,1:3);
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
%line(xu,yu,'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor',grey,'MarkerEdgeColor',grey), hold on

ind=find(xu < maxlon & xu > minlon & yu < maxlat & yu > minlat );

xlin=[xu yu]; 

D = 0.3; % distance threshold for selecting nearby points
style = ':';
color = 'b';
width = 2;
color = [rand(1) rand(1) rand(1)];
X=[]; Y=[]; Len=[]; S=[];

f = 'fc2.new';
f = 'fc2.dat';
%f = 'lin_segments.dat';
%f = 'lin_segments.rev';
if isfile(f)
 fid = fopen(f,'r');
 load(f);
% x=lin_segments;
 x=fc2;
 [m,n]=size(x);
 e=x(:,n);  % errors
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

% for i=1:m
 % line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
% end
 fclose(fid);
end

%pause

color = [0 0 0];  
fid=fopen('junk','wt');

figure('Renderer', 'painters', 'Position', [10 10 600 900])
clf

% % plot all of the clusters:
x0=X; y0=Y; L=Len; strike=S;
m=5; n=5;
m=4; n=5;
m=14; n=5;

%ha = tight_subplot(m,n,[-0.17 .045],[0 -0.08],[.04 .01]); % 5x5
%ha = tight_subplot(m,n,[0.01 .045],[0.62 0.01],[.04 .01]);  % 2x5 
ha = tight_subplot(m,n,[0.01 .045],[0.25 0.01],[.04 .01]);  % 3x5 

%for i=1:n*m  % use for printing panels

for i=1:70 % do all clusters
%for i=1:65
  j=i;
%  j=i+n*m;
%  j=i+50;
%for i=1:length(x0)
%for i=length(x0):-1:1

% sp(i)=subplot(m,n,i);
 axes(ha(i)); 

 %figure(1)
 % figure('Renderer', 'painters', 'Position', [10 10 200 200])
 % clf

 [xmin,xmax,ymin,ymax]=plot_fm(i,x0(j),y0(j),xlin,FM,strike(j),L(j),D,color,style,width,fid,qn,e);
 %tit=[cellstr(['segment number : ' num2str(i)])];
 %title(char(tit),'Fontsize',14);
 axis equal;
 set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax],'box','on','FontSize',Fs);
 set(gca,'box','on','FontSize',10);
% saveas(gcf, ['streak_fm_' num2str(i)], 'png')
 %pause
end

% set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

% saveas(gcf, ['streak_fm'], 'png')
print(gcf,'-dpng','-r300','streak_fm1');

fclose(fid);
return


function [xmin,xmax,ymin,ymax]=plot_fm(k,x0,y0,xlin,FM,strike,L,D,color,style,width,fid,qn,e)

%[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,fid);
[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,'r','-',3,0);
%plot([x1 x2],[y1 y2],'ko'), hold on
v1=[x1 y1 0];  v2=[x2 y2 0]; 
P(1)=(y2-y1)/(x2-x1);
P(2)=y1-P(1)*x1;

nx=[]; ny=[];
pts=[xlin(:,1) xlin(:,2) 0*xlin(:,1)];

%[ind] = point_to_line(pts, v1, v2, D);
if D < L/3
 distp=D;
else
 distp=L/3;
end

[ind] = point_to_line(pts, v1, v2, distp);
nx = xlin(ind,1);
ny = xlin(ind,2);
Q = mean(qn(ind));

xmin=min(nx); xmax=max(nx); 
ymin=min(ny); ymax=max(ny);
dim=max(xmax-xmin,ymax-ymin);
%dim2=min(xmax-xmin,ymax-ymin);
mar=2;  % factor for making margins
xmin=xmin-dim/mar;  xmax=xmax+dim/mar;
ymin=ymin-dim/mar;  ymax=ymax+dim/mar;

 if (k==1 | k==6) % custom bounds to better position certain sub-plots
%  ymin=ymin+dim/mar/2;  ymax=ymax-dim/mar/2;
 end

dim2=sqrt((xmax-xmin)^2+(ymax-ymin)^2); % factor for scaling the beach balls
  
% all events
line(xlin(:,1),xlin(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor',[0.6 0.6 0.6]), hold on

% cluster events
line(nx,ny,'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k'), hold on

[Err]=err_2D_data(k, nx,ny, P);
%Err=e(k);

if k==2
%  Err=8;
end

%[x1,y1,x2,y2]=plot_line(x0,y0,strike+Err,L,'m','-',1,0);
%[x1,y1,x2,y2]=plot_line(x0,y0,strike-Err,L,'b','-',1,0);

fm=FM(ind,1:3);

arr=[nx ny fm];

deg2rad=pi/180.;
comp=zeros(1,6);
sm=[];
num_events=length(arr(:,1));
 for i=1:num_events
  [fm]=sdr2mij(arr(i,3:5),deg2rad);
  comp=comp+fm;
  % this is to estimate variability:
  [s1,d1,s2,d2] = mij2sd(fm);
  sm=[sm min([s1 s2])];
 end
stdm=std(sm); 

%fid=fopen('junk1.dat','wt');
  if fid > 0
%   fprintf(fid,'%8.4f %8.4f %8.4f %8.4f %f %f %f %f %f %f %f\n',x1,y1,x2,y2,comp(:),Err);
   fprintf(fid,'%f\n',Err);
  end
%fclose(fid);
 
[s1,d1,s2,d2] = mij2sd(comp);

% % average mechanism:
style = '--';
%color = 'b';
width = 3;
x0=mean(arr(:,1));
y0=mean(arr(:,2));
%ds=rem(abs(strike-s1),90);
ds=abs(strike-s1);
if abs(rem(round(ds/90),2))==0
%if ds<45
  s=s1;
else
  s=s2;
end
%[x1,y1,x2,y2]=plot_line(x0,y0,s,L,color,style,width,fid);
%[x1,y1,x2,y2]=plot_line(x0,y0,s,L,color,style,width,0);

del=dim/2;

% scale (1 km):
%line([xmin+0.3 xmin+1.3],[ymax-0.2 ymax-0.2],'Color','m','LineStyle','-','LineWidth',3);

%line([mean(nx) mean(nx)+del],[mean(ny) mean(ny)-del],'Color','k','LineStyle','-','LineWidth',0.5);
%bb(comp, mean(nx)+0.9*del, mean(ny)-0.3*del, dim/5, 0, 'b'), hold on
focalmech(comp, mean(nx)+0.9*del, mean(ny)-0.3*del, dim2/6, 'b','dc'), hold on
text(mean(nx)+0.7*del, mean(ny)-1.1*del,num2str(num_events),'FontSize',12), hold on
text(mean(nx)-0.1*del, mean(ny)+1.3*del,num2str(Err,'%4.0f\n'),'FontSize',14,'Color','r'), hold on
text(xmin+del/3, ymax-del/3, ['(' char(k + 96) ')'],'FontSize',12), hold on
%text(xmin+del/10, ymax-del/10, num2str(k),'FontSize',14), hold on
%text(mean(nx)+del, mean(ny)-1.8*del,num2str(k),'FontSize',12), hold on
% text(mean(nx)+0.5, mean(ny)+0.5,num2str(k),'FontSize',11), hold on
% text(mean(nx)+0.57, mean(ny)+0.57,num2str(Q),'FontSize',11,'Color',color), hold on
% text(mean(nx)+1.2, mean(ny)+0.59,num2str(stdm),'FontSize',14,'Color','k'), hold on
% text(mean(nx)+1.2, mean(ny)+0.39,num2str(mean(sm)),'FontSize',14,'Color','k'), hold on
% text(mean(nx)+1.2, mean(ny)+0.19,num2str(s1),'FontSize',14,'Color','k'), hold on
%text(mean(nx)+0.5, mean(ny)-0.5,num2str(strike),'FontSize',11,'Color',color), hold on
%text(mean(nx)+0.5, mean(ny)-0.53,num2str(s1),'FontSize',11,'Color',color), hold on
%text(mean(nx)+0.5, mean(ny)-0.56,num2str(ds),'FontSize',11,'Color',color), hold on

%set(gca,'xlim',[xmin xmax],'ylim',[ymin ymax],'box','on','FontSize',12);

end

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


  line([x1 x2],[y1 y2],'lineWidth',width,'lineStyle',style,'Color',color), hold on
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


function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
%function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
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
