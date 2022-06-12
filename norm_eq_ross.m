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

figure(1)

load eqn1.ross; % iteratively selected aftershocks (remove outliers)
x=eqn1;
%[x1,y1]=utm2ll(x(:,1),x(:,2),0,1);
%x1=(x1-xo)*1e-3;
%y1=(y1-yo)*1e-3;

line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor','k','MarkerEdgeColor','k'), hold on

load eqn2.ross; % iteratively selected aftershocks (remove outliers)
x=eqn2;
%line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor','r','MarkerEdgeColor','r'), hold on

set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat]);
%return


[X,Y]=utm2ll(x(:,1),x(:,2),0,1);
%[X,Y]=utm2ll(x,y,0,1);
%X=x(:,1);
%Y=x(:,2);
X=(X-xo)*1e-3;
Y=(Y-yo)*1e-3;
xl=[min(X) max(X)];
yl=[min(Y) max(Y)];
xval=xl(1):(xl(2)-xl(1))/100:xl(2);

X=X/(xl(2)-xl(1));
Y=Y/(xl(2)-xl(1));
X=X-min(X);
Y=Y-min(Y);
[x1,ind]=sort(X);
y1=Y(ind);
new = [x1,y1];

%save clust3.nrm new -ascii
%save lin.utm new -ascii
save eqn2.nrm new -ascii

return

X=x(:,1);
Y=x(:,2);
event = find(X> minlon & X < maxlon & Y < maxlat & Y > minlat);
x=X(event); 
y=Y(event);
new = [x,y];

save lin.ll new -ascii
