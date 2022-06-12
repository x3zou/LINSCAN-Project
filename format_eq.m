clear all;
close all;
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

% regional
minlon = -118.15; maxlon = -117.2; minlat = 35.45; maxlat = 36.2; 

% clust1
%minlon = -118.15; maxlon = -117.95; minlat = 35.6; maxlat = 35.85; 

% clust2
%minlon = -117.64; maxlon = -117.603; minlat = 35.65; maxlat = 35.69; 

% clust3
%minlon = -117.6; maxlon = -117.52; minlat = 35.607; maxlat = 35.66; 

% clust4
%minlon = -117.6; maxlon = -117.52; minlat = 35.78; maxlat = 35.82; 


%load Ridgecrest1.trace;
%x=Ridgecrest1;
%load eq_lin.dat;
%x=eq_lin;
load eq1.dat; % aftershocks of M7.1
x=eq1;

%load eq1.dat; % aftershocks of M7.1
%x=eq1;
%load m6.aft; % aftershocks of M6.4
%x=m6;

figure(1)
line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','k'), hold on

load m6.aft; % aftershocks of M7.1
x=m6;
%load Ridgecrest_fault_ll.dat;
%x=Ridgecrest_fault_ll;
line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','r'), hold on

set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat]);
return

X=x(:,1);
Y=x(:,2);
event = find(X> minlon & X < maxlon & Y < maxlat & Y > minlat);
x=X(event); 
y=Y(event);
new = [x,y];

save lin.ll new -ascii

%[X,Y]=utm2ll(x(:,1),x(:,2),0,1);
[X,Y]=utm2ll(x,y,0,1);
%X=x(:,1);
%Y=x(:,2);
X=(X-xo)*1e-3;
Y=(Y-yo)*1e-3;
xl=[min(X) max(X)];
yl=[min(Y) max(Y)];
xval=xl(1):(xl(2)-xl(1))/100:xl(2);

%X=X/(xl(2)-xl(1));
%Y=Y/(xl(2)-xl(1));
%X=X-min(X);
%Y=Y-min(Y);
[x1,ind]=sort(X);
y1=Y(ind);
new = [x1,y1];

%save clust3.nrm new -ascii
save lin.utm new -ascii

return

unix('cp eq1.dat eqn1.dat');
NumIter=4; %number of iterations


for j=1:NumIter

load eqn1.dat; % iteratively selected aftershocks (remove outliers)
x=eqn1;
[x1,y1]=utm2ll(x(:,1),x(:,2),0,1);

x1=(x1-xo)*1e-3;
y1=(y1-yo)*1e-3;

[xs,I]=sort(x1);
ys=y1(I);

% remove non-increasing data points (twice):
s=find(diff(xs)==0);
L = ~ismember(I,s);
xs=x1(I(L));
ys=y1(I(L));
s=find(diff(xs)==0);
length(s)
I=find(diff(xs)~=0);
xs=xs(I);
ys=ys(I);
s=find(diff(xs)==0);
length(s)
I=find(diff(xs)~=0); 

win = 1; % window size, km 
TF = isoutlier(ys,'movmedian',win,'SamplePoints',xs);
length(find(TF~=0))
plot(xs(TF),ys(TF),'xr'), hold on
%plot(xs(~TF),ys(~TF),'.b'), hold on
%legend('Data','Outlier')

% convert back to lat/lon
UTM=11; % UTM zone
[x1,y1]=utm2ll(xs(~TF)*1e3+xo,ys(~TF)*1e3+yo,UTM,2);

new = [x1,y1];

save eqn1.dat new -ascii

end %iterations

%return

x1=xs(~TF);
y1=ys(~TF);
new = [x1,y1,0*x1];

% just to get points on a line:
[E,P]=fit_2D_data(x1,y1,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;

v1=[xval(1),yval(1),0]; 
v2=[xval(length(xval)),yval(length(yval)),0]; 
perp = dist_point_to_line(new,v1,v2);

dist = 5; % distance threshold, km
s=find(perp < dist);
x=x1(s); 
y=y1(s); 

plot(X,Y,'.g'), hold on
plot(x,y,'.m'), hold on

new = [x,y];

save eqn2.dat new -ascii

%return

% fit straight lines:

[E,P]=fit_2D_data(x,y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','k','LineWidth',2), hold on

[p,s]=polyfit(x,y,1);
yval=p(2)+p(1)*xval;
line(xval,yval,'LineStyle','--','Color','k','LineWidth',1), hold on

[pj,s]=polyfit(y,x,1);
p2(1)=1/pj(1);
p2(2)=-pj(2)/pj(1);
yval=p2(2)+p2(1)*xval;
line(xval,yval,'LineStyle',':','Color','k','LineWidth',1), hold on

[b,bintr,bintjm] = gmregress(x,y);
A0=b(1); S0=b(2);
%[s0]=lin_fit(x,y);
%A0=mean(Y(:))+1/s0*mean(X(:)); % intercept
%A0=mean(y(:))-s0*mean(x(:)); % intercept
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','c','LineWidth',2), hold on

x=X; y=Y;
[E,P]=fit_2D_data(x,y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','b','LineWidth',2), hold on

[p,s]=polyfit(x,y,1);
yval=p(2)+p(1)*xval;
line(xval,yval,'LineStyle','--','Color','b','LineWidth',1), hold on

[pj,s]=polyfit(y,x,1);
p2(1)=1/pj(1);
p2(2)=-pj(2)/pj(1);
yval=p2(2)+p2(1)*xval;
line(xval,yval,'LineStyle',':','Color','b','LineWidth',1), hold on

load Ridgecrest1.trace;  %M7
x=Ridgecrest1;
[X,Y]=utm2ll(x(:,1),x(:,2),0,1);
X=(X-xo)*1e-3;
Y=(Y-yo)*1e-3;
plot(X,Y,'oy'), hold on

[E,P]=fit_2D_data(X,Y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','r','LineWidth',2), hold on
atan(s0)/pi*180

axis('equal');

figure(2)

load eq2.dat; % aftershocks of M6.4
x=eq1;
[X,Y]=utm2ll(x(:,1),x(:,2),0,1);
%X=x(:,1);
%Y=x(:,2);
X=(X-xo)*1e-3;
Y=(Y-yo)*1e-3;
xl=[min(X) max(X)];
yl=[min(Y) max(Y)];
xval=xl(1):(xl(2)-xl(1))/100:xl(2);

unix('cp m6.aft eqf1.dat');
%return

NumIter=4; %number of iterations


for j=1:NumIter

load eqf1.dat; % iteratively selected aftershocks (remove outliers)
x=eqf1;
[x1,y1]=utm2ll(x(:,1),x(:,2),0,1);

x1=(x1-xo)*1e-3;
y1=(y1-yo)*1e-3;

[xs,I]=sort(x1);
ys=y1(I);

% remove non-increasing data points (twice):
s=find(diff(xs)==0);
L = ~ismember(I,s);
xs=x1(I(L));
ys=y1(I(L));
s=find(diff(xs)==0);
length(s)
I=find(diff(xs)~=0);
xs=xs(I);
ys=ys(I);
s=find(diff(xs)==0);
length(s)
I=find(diff(xs)~=0); 

win = 1; % window size, km 
TF = isoutlier(ys,'movmedian',win,'SamplePoints',xs);
length(find(TF~=0))
plot(xs(TF),ys(TF),'xr'), hold on
%plot(xs(~TF),ys(~TF),'.b'), hold on
%legend('Data','Outlier')

% convert back to lat/lon
UTM=11; % UTM zone
[x1,y1]=utm2ll(xs(~TF)*1e3+xo,ys(~TF)*1e3+yo,UTM,2);

new = [x1,y1];

save eqf1.dat new -ascii

end %iterations

%return

x1=xs(~TF);
y1=ys(~TF);
new = [x1,y1,0*x1];

% just to get points on a line:
[E,P]=fit_2D_data(x1,y1,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;

v1=[xval(1),yval(1),0]; 
v2=[xval(length(xval)),yval(length(yval)),0]; 
perp = dist_point_to_line(new,v1,v2);

dist = 5; % distance threshold, km
s=find(perp < dist);
x=x1(s); 
y=y1(s); 

plot(X,Y,'.g'), hold on
plot(x,y,'.m'), hold on

new = [x,y];

save eqf2.dat new -ascii

%return

% fit straight lines:

[E,P]=fit_2D_data(x,y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','k','LineWidth',2), hold on

[p,s]=polyfit(x,y,1);
yval=p(2)+p(1)*xval;
line(xval,yval,'LineStyle','--','Color','k','LineWidth',1), hold on

[pj,s]=polyfit(y,x,1);
p2(1)=1/pj(1);
p2(2)=-pj(2)/pj(1);
yval=p2(2)+p2(1)*xval;
line(xval,yval,'LineStyle',':','Color','k','LineWidth',1), hold on

[b,bintr,bintjm] = gmregress(x,y);
A0=b(1); S0=b(2);
%[s0]=lin_fit(x,y);
%A0=mean(Y(:))+1/s0*mean(X(:)); % intercept
%A0=mean(y(:))-s0*mean(x(:)); % intercept
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','c','LineWidth',2), hold on

x=X; y=Y;
[E,P]=fit_2D_data(x,y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','b','LineWidth',2), hold on

[p,s]=polyfit(x,y,1);
yval=p(2)+p(1)*xval;
line(xval,yval,'LineStyle','--','Color','b','LineWidth',1), hold on

[pj,s]=polyfit(y,x,1);
p2(1)=1/pj(1);
p2(2)=-pj(2)/pj(1);
yval=p2(2)+p2(1)*xval;
line(xval,yval,'LineStyle',':','Color','b','LineWidth',1), hold on


load Ridgecrest2.trace;  %M6
x=Ridgecrest2;
[X,Y]=utm2ll(x(:,1),x(:,2),0,1);
X=(X-xo)*1e-3;
Y=(Y-yo)*1e-3;
plot(X,Y,'ob'), hold on

[E,P]=fit_2D_data(X,Y,'no'); % total least squares
s0=P(1); A0=P(2);
yval=A0+s0*xval;
line(xval,yval,'LineStyle','-','Color','r','LineWidth',2), hold on
atan(s0)/pi*180


%xo=orig(1); yo=orig(2);
%[x1,y1]=utm2ll(xo,yo,0,1);
%xo=orig(1); yo=orig(2)+0.1;
%[x2,y2]=utm2ll(xo,yo,0,1);
%len = sqrt((x2-x1)^2+(y2-y1)^2)
%a = asin((x2-x1)/len)/pi*180

axis('equal');

return


function perp = dist_point_to_line(pt, v1, v2)
 % Compute distance
 v1 = repmat(v1,size(pt,1),1);
 v2 = repmat(v2,size(pt,1),1);
 a = v1 - v2;
 b = pt - v2;
 perp = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));
 %a = v1 - v2;
 %b = pt - v2;
 %perp = norm(cross(a,b)) / norm(a);
end
