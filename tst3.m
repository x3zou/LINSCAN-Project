clear all;
close all;

x=-10:0.1:10;
a=-10.3; b=0.7;
y=a*x+b;

R=20; % scale for scatter
xh=R*rand(size(x))+x;
yh=R*rand(size(y))+y;

plot(xh,yh,'k.'), hold on
plot(x,y,'r-'), hold on

[x1,y1,x2,y2]=lin_fit2(xh,yh);
line([x1 x2],[y1 y2], 'lineStyle','-'), hold on

return

x1=min(x); x2=max(x); 
y1=min(y); y2=max(y); 
P(1)=(y2-y1)/(x2-x1);
P(2)=y1-P(1)*x1;

%[Err]=err_2D_data(xh,yh, P)

mdlr = fitlm(xh,yh,'RobustOpts','on');
mdl = fitlm(xh,yh);
%plotResiduals(mdl,'probability')


beta1 = mdl.Coefficients.Estimate;
y1=xh*beta1(2)+beta1(1);
plot(xh,y1,'g*'), hold on
  beta2 = mdlr.Coefficients.Estimate;
y2=xh*beta2(2)+beta2(1);
plot(xh,y2,'k--'), hold on
err2 = mdlr.Coefficients.SE

y3=xh*(beta2(2)+err2(2))+beta2(1);
plot(xh,y3,'r--'), hold on
y4=xh*(beta2(2)-err2(2))+beta2(1);
plot(xh,y4,'b--'), hold on


% try later for fitting a 3D plane:
%  X=[xh',y2'];
%  b = robustfit(X,yh');
%y3=xh*b(2)+b(1);

%plot(xh,y3,'bs'), hold on

