clear all;
close all;

minlon = -60; maxlon = 25; minlat = 10; maxlat = 75; 
% %minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
% %minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
Fs=14;

pr=1;
fault_tr = 1; % flag for plotting the fault trace
% regional
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
 grey = [1 1 1]*0.5;
% %grey = [0.4,0.4,0.4];
% 
%load lin.utm;
%X=lin;
load fm.utm;
X=fm;


ind=find(X(:,2)<10);
X(ind,:)=[];
% 
% figure(1)
% %line(xlin(:,1),xlin(:,2),'lineStyle','none','Marker','.','MarkerSize',0.1,'MarkerFaceColor',grey,'MarkerEdgeColor',grey), hold on
% 
scatter(X(:,1),X(:,2),1,grey,'.'), hold on
% 
f = 'lin_segments.dat';
if isfile(f)
 fid = fopen(f,'r');
 load lin_segments.dat;
 [m,n]=size(lin_segments);
 for i=1:m
  line([lin_segments(i,1) lin_segments(i,3)],[lin_segments(i,2) lin_segments(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
 end
 fclose(fid);
end
%return

pause
% 
% % DBSCAN:
% 
minpts = 15; % Minimum number of neighbors for a core point

% kD = pdist2(X,X,'euc','Smallest',minpts); % The minpts smallest pairwise distances
% 
% figure(2)
% plot(sort(kD(end,:)));
% title('k-distance graph')
% %title('k-distance graph %d',minpts)
% xlabel('Points sorted with 50th nearest distances')
% ylabel('N-th nearest distances')
% grid
% 
f = 'lin_segments.dat';
fid = fopen(f,'wt');

cnt=0;
de=[0.5 0.5 0.7 1 1.5]; % increments for epsilon 
amblims = [min(X(:,1)) max(X(:,1)); min(X(:,2)) max(X(:,2))];

while cnt < 5
    cnt=cnt+1;
    
epsilon = 0.0 + de(cnt);  % distance that defines neighbors in a cluster
%minpts=minpts-1;
% 
labels = dbscan(X,epsilon,minpts);

%clusterer = clusterDBSCAN('MinNumPoints',minpts,'Epsilon',epsilon);
%labels=clusterer(X);
%clusterer = clusterDBSCAN('MinNumPoints',minpts,'Epsilon',epsilon, ...
%'EnableDisambiguation',true,'AmbiguousDimension',[1 2]);
%labels=clusterer(X,amblims);
% 

figure(3)
clf
gscatter(X(:,1),X(:,2),labels);
hLeg= findobj(gcf,'tag','legend'); 
set(hLeg,'visible','off');
%set(l,'location','northeastoutside');
%title('epsilon = 2 and minpts = 10')
tit=[cellstr(['number of points: ' num2str(length(X(:,1))) ])];
title(char(tit),'Fontsize',14);
grid
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
if pr==1
 print(gcf, '-djpeg90', '-r300', ['clust1' char(num2str(cnt))])
end
%pause

figure(5)
clf
ind=find(labels>0);
%  scatter(X(ind,1),X(ind,2),1,'r','.');
gscatter(X(ind,1),X(ind,2),labels(ind));
hLeg= findobj(gcf,'tag','legend'); 
set(hLeg,'visible','off');
tit=[cellstr(['number of clusters: ' num2str(max(labels)) ])];
title(char(tit),'Fontsize',14);
grid
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
if pr==1
 print(gcf, '-djpeg90', '-r300', ['clust2' char(num2str(cnt))])
end
%pause

% now run RANSAC to get linear trends out of identified clusters:
% RANSAC:

sampleSize = 2; % number of points to sample per trial
maxDistance = 2; % max allowable distance for inliers
mclust=[]; % combined indecies of selected points

denom=12;
min_size=15;

%figure(6)
clf
for i=1:max(labels)
%%for i=1:0
 points=[];
 clust=find(labels==i);
  if length(clust)>min_size
   points = [X(clust,1) X(clust,2)];

%   fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
%   evalLineFcn = ...   % distance evaluation function
%   @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);
% 
%   [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
%                                     sampleSize,maxDistance,'Confidence',90);
%   modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);
% 
%   inlierPts = points(inlierIdx,:);
%   x = [min(inlierPts(:,1)) max(inlierPts(:,1))];
%   y = modelInliers(1)*x + modelInliers(2);
%   v1=[x(1) y(1) 0];
%   v2=[x(2) y(2) 0];
  
   scatter(points(:,1),points(:,2),3,'m','o'), hold on
   [Err, P, C] = fit_2D_data(points(:,1),points(:,2), 'no');
   x=[C(1) C(3)]; 
   y=[C(2) C(4)];
   v1=[x(1) y(1) 0];
   v2=[x(2) y(2) 0];
   %plot(x, y, 'c-','LineWidth',3)
  
  %pause
%  [ind2,dist] = point_to_line([points(:,1) points(:,2)], v1, v2, 5*maxDistance);
   [ind2,dist] = point_to_line([points 0*points(:,1)], v1, v2, 5*maxDistance);
 % fprintf('dist: %f, length: %f \n',mean(dist),norm(v1 - v2)/10);
%   if mean(dist)>norm(v1 - v2)/(18); % condition for accepting a good linear fit
   if mean(dist)>norm(v1 - v2)/(denom+cnt); % condition for accepting a good linear fit
%  if median(dist)>norm(v1 - v2)/(denom-cnt); % condition for accepting a good linear fit
%  if max(dist)>norm(v1 - v2)/(5-0.2*cnt); % condition for accepting a good linear fit
    plot(x, y, 'g-','LineWidth',3)
%   text(mean(x)+0.3, mean(y)-0.3, char(num2str(i)),'Color','g')
%    fprintf(fid,'%9.4f %9.4f %9.4f %d %s\n',max(dist)/norm(v1 - v2), median(dist)/norm(v1 - v2), mean(dist)/norm(v1 - v2),i,'f');
   else
    plot(x, y, 'k-','LineWidth',3)
%   text(mean(x)+0.3, mean(y)-0.3, char(num2str(i)),'Color','k')
    fprintf(fid,'%9.4f %9.4f %9.4f %9.4f\n',x(1), y(1), x(2), y(2));
%   fprintf(fid,'%9.4f %9.4f %9.4f %d %s\n',max(dist)/norm(v1 - v2), median(dist)/norm(v1 - v2), mean(dist)/norm(v1 - v2),i,'p');
   % flag accepted points for removal from the catalog for subsequent iterations:
    mclust=[mclust' clust']';
   end % if mean(dist)>norm(v1 - v2)/(denom-cnt)
  %legend('Noisy points','Least squares fit','Robust fit');
  %hold off
  end % if length(clust)>min_size
 end % of loop over identified clusters
  X(mclust,:)=[];
%pause
 if pr==1
  grid
  set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
  print(gcf, '-djpeg90', '-r300', ['clust3' char(num2str(cnt))])
 end
end % of iterations

fclose(fid);

% plot the selected fits:
X=fm;
ind=find(X(:,2)<10);
X(ind,:)=[];
% 
figure(7)
scatter(X(:,1),X(:,2),1,grey,'.'), hold on
% 
f = 'lin_segments.dat';
if isfile(f)
 fid = fopen(f,'r');
 load lin_segments.dat;
 [m,n]=size(lin_segments);
 for i=1:m
  line([lin_segments(i,1) lin_segments(i,3)],[lin_segments(i,2) lin_segments(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
 end
 fclose(fid);
end

return
% 
% return

%f = 'sc2019_hash_ABCD_so.focmec.scedc.txt';
f = 'sc.focmec.scedc.txt';
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
%Q=L(ind);
%M=[x y];
%save jjj M -ascii;

%return

FM=FM(ind,1:3);

[xx,yy]=utm2ll(x(:),y(:),0,1);
% !!!!!
xshift=0.0; % km; stress drop locations are systematically to the E from Lin relocated catalog

xu=reshape((xx-xo)*1e-3,size(x)) - xshift;
yu=reshape((yy-yo)*1e-3,size(x));
line(xu,yu,'lineStyle','none','Marker','.','MarkerSize',3,'MarkerFaceColor','m','MarkerEdgeColor','m'), hold on

ind=find(xu < maxlon & xu > minlon & yu < maxlat & yu > minlat );
xlin=[xu(ind) yu(ind) z(ind)]; 
save jjj xlin -ascii;
return

xlin=[xu yu]; 

strike = 42;
style = ':';
color = 'b';
width = 2;

%fid=fopen('fc.dat','wt');

D = 0.7; % distance threshold for selecting nearby points

color = [rand(1) rand(1) rand(1)];

X=[]; Y=[]; Len=[]; S=[];

% % Left-lateral:
% fid=fopen('LL.dat','wt');
% 
% %x0 =-16.6; y0 = 63; L = 3;  
% x0 =-17.9; y0 = 62; L = 5;  
% 
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %plot_fm(x0,y0,xlin,FM,strike,L,D,color,style,width,fid);
% 
% % Do a bunch:
% 
% 
% x0 =-50; y0 = 25; L = 20;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% %plot_line(x0,y0,strike,L,color,style,width,fid);
% 
% x0 =-11; y0 = 19; L = 9;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-4; y0 = 35; L = 7;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-7; y0 = 60; L = 9;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-3.5; y0 = 58; L = 9;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-12.751225; y0 = 62.837010; L = 17;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %x0 =-36.158088; y0 = 69.699755; L = 5;  % normal mechanism
% x0 =-34.2; y0 = 67.7; L = 5;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %x0 =-31; y0 = 63; L = 15;
% %X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %x0 =-31; y0 = 63; L = 15;
% %X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =17; y0 = 52.5; L = 6;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-5.5; y0 = 15; L = 7;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% % plot all of the clusters:
% x0=X; y0=Y; L=Len; strike=S;
% for i=1:length(x0)
%  color = [rand(1) rand(1) rand(1)];   
%  plot_fm(x0(i),y0(i),xlin,FM,strike(i),L(i),D,color,style,width,fid);
% end
% 
% fclose(fid);
% %return
% 
% 
%  fid=fopen('RL.dat','wt');
%  X=[]; Y=[]; Len=[]; S=[];
% 
% % Right-lateral:
% 
% strike = -40;
% %style = ':'; color = 'r'; color = 'k';
% 
% 
% x0 =-20.5; y0 = 29.2; L = 3; 
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0=-40; y0= 48.28; L = 10;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %x0=-12.230392; y0= 27.787990; L = 6;  % mostly normal
% %X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0=10.379902; y0=53.768382; L = 7;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-14; y0 = 51; L = 10;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-25.5; y0 = 55; L = 10;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =-2; y0 = 60; L = 5;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% x0 =19; y0 = 50; L = 5;
% X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% %x0 =-47; y0 = 31; L = 5;
% %x0 =-12; y0 = 31.8; L = 3;  % mostly normal
% %X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 
% 
% % plot all of the clusters:
% x0=X; y0=Y; L=Len; strike=S;
% for i=1:length(x0)
%  color = [rand(1) rand(1) rand(1)];   
%  plot_fm(x0(i),y0(i),xlin,FM,strike(i),L(i),D,color,style,width,fid);
% end
% 
% fclose(fid);
% 
% return
% 

fid=fopen('other.dat','wt');

strike = 35;

x0 =-14.85; y0 = 69; L = 2; 
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

strike = 24;

x0 =-17.1; y0 = 67; L = 4; 
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

strike = 32;

x0 =-7.8; y0 = 62.4; L = 2; 
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

strike = -15;
x0 =-32.8; y0 = 53.7; L = 2; 
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

strike=-23;
x0 =19; y0 = 52.9; L = 4; % match Lin locations
x0 =19.7; y0 = 52.9; L = 4; % match stress drop catalog locations (0.5 km shift to E)
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

strike=-23;
x0 =-8.2; y0 = 19; L = 2; 
X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

%strike=-20;
%x0 =-13.3; y0 = 65; L = 2;
%X=[X x0]; Y=[Y y0]; Len=[Len L]; S=[S strike]; 

% plot all of the clusters:
x0=X; y0=Y; L=Len; strike=S;
for i=1:length(x0)
 color = [rand(1) rand(1) rand(1)];   
 plot_fm(x0(i),y0(i),xlin,FM,strike(i),L(i),D,color,style,width,fid);
end

fclose(fid);

if fault_tr==1

 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');

end

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

% f = 'fm.pcl';
% savepcd(f,X');
% f = '/Users/fialko/Downloads/ji/pointclouds/room.con';
% [y]=loadpcd(f);
% size(y)
% 
% return
% 
% ptCloud = pointCloud(X);
% figure
% pcshow(ptCloud);
% maxDistance = 1;
% 
% roi = [-inf,inf;-inf,inf;-inf,inf];
% sampleIndices = findPointsInROI(ptCloud,roi);
% ptCloudSegm = select(ptCloud,sampleIndices,'OutputSize','full');
% distThreshold = 0.3;
% [labels,numClusters] = pcsegdist(ptCloudSegm,distThreshold);
% labelColorIndex = labels+1;
% pcshow(ptCloud.Location,labelColorIndex)
% colormap([hsv(numClusters);[0 0 0]])
% title('Point Cloud Clusters');
% return
% 
% roi = [-60,-40;10,40;-15,0];
% sampleIndices = findPointsInROI(ptCloud,roi);
% [model2,inlierIndices,outlierIndices] = pcfitplane(ptCloud,...
%             maxDistance,'SampleIndices',sampleIndices);
% plane2 = select(ptCloud,inlierIndices);
% figure
% pcshow(plane2)
% title('First Plane')
% return

% % this is input for sdr2mij.py:
% fid=fopen('sdr.dat','wt');
% for i=1:length(arr(:,1))
%  fprintf(fid,'%f %f %d %d %d\n',arr(i,1:2),arr(i,3:5));
% end
% fclose(fid);
% 
% 
% % this is output of sdr2mij.py:
% load mij.dat;

%[fm]=sdr2mij(s,d,r,pi/180.)

% nx=mij(:,1); 
% ny=mij(:,2);
% 
% fm= mij(:,3:8);
% %return
% X=nx*1e3+xo;
% Y=ny*1e3+yo;
% [nx,ny]=utm2ll(X,Y,11,2);
% 
% % figure out the correct lat/lon aspect ratio
% del=0.01;
% orx=mean(nx(:)); ory=mean(ny(:));
% [tx1,ty1]=utm2ll(orx,ory,0,1);
% [tx2,ty2]=utm2ll(orx+del,ory,0,1);
% [tx3,ty3]=utm2ll(orx,ory+del,0,1);
% llrat=[1 abs((tx2-tx1)/(ty3-ty1)) 1];

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

function plot_fm(x0,y0,xlin,FM,strike,L,D,color,style,width,fid)

[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,fid);
plot([x1 x2],[y1 y2],'ko'), hold on
v1=[x1 y1 0];  v2=[x2 y2 0]; 

nx=[]; ny=[];
pts=[xlin(:,1) xlin(:,2) 0*xlin(:,1)];

[ind] = point_to_line(pts, v1, v2, D);
nx = xlin(ind,1);
ny = xlin(ind,2);
line(nx,ny,'lineStyle','none','Marker','.','MarkerSize',7,'MarkerFaceColor',color,'MarkerEdgeColor',color), hold on

[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,'k',style,4,fid);

fm=FM(ind,1:3);

arr=[nx ny fm];

deg2rad=pi/180.;
comp=zeros(1,6);
num_events=length(arr(:,1));
 for i=1:num_events
  [fm]=sdr2mij(arr(i,3:5),deg2rad);
  comp=comp+fm;
 end

 
[s1,d1,s2,d2] = mij2sd(comp);

% % average mechanism:
style = '--';
%color = 'b';
width = 3;
x0=mean(arr(:,1));
y0=mean(arr(:,2));
%[x1,y1,x2,y2]=plot_line(x0,y0,s2,L,color,style,width,fid);

bb(comp, mean(nx), mean(ny), 0.7, 0, color), hold on
text(mean(nx)+1.5, mean(ny)-1.5,num2str(num_events),'FontSize',14), hold on

end
