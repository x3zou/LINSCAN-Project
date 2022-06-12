close all;
clear all;

y=[0:100];
x=0*y;
x=[0:100];
y=0*x;
A=10;
del=A*rand(size(y));
Y=y+del;

del=A*rand(size(y));
X=x+del;

figure(1)
plot(X,Y,'k.'), hold on

axis equal;

[Err, P, C] = fit_2D_data(X, Y, 'no');
axis equal;
line([C(1) C(3)], [C(2) C(4)], 'Color','g','LineWidth',3), hold on 

%set(gca,'xlim',[min(X) max(X)],'ylim',[min(Y) max(Y)]);
%return

%del2=(max(X)-min(X))/10;
%x1=min(X):del2:max(X);
%y1=P(2)*x1+P(1);
%plot(x1,y1,'r*'), hold on
%y2=P(2)*X+P(1);
%plot(X,y2,'ms'), hold on
%axis equal;
%del2=del2*3;
%set(gca,'xlim',[min(X)-del2 max(X)+del2],'ylim',[min(Y)-del2 max(Y)+del2]);

%return


sampleSize = 2; % number of points to sample per trial
maxDistance = 2; % max allowable distance for inliers
  points = [X' Y'];
  fitLineFcn = @(points) polyfit(points(:,1),points(:,2),1); % fit function using polyfit
  evalLineFcn = ...   % distance evaluation function
  @(model, points) sum((points(:, 2) - polyval(model, points(:,1))).^2,2);

  [modelRANSAC, inlierIdx] = ransac(points,fitLineFcn,evalLineFcn, ...
                                    sampleSize,maxDistance,'Confidence',90);
  modelInliers = polyfit(points(inlierIdx,1),points(inlierIdx,2),1);

  inlierPts = points(inlierIdx,:);
  x3 = [min(inlierPts(:,1)) max(inlierPts(:,1))];
  y3 = modelInliers(1)*x3 + modelInliers(2);

plot(x3,y3,'b-'), hold on

p = polyfit(X,Y,1);
f = polyval(p,X);

plot(X,f,'r-'), hold on

