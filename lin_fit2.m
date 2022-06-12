%function [x1,y1,x2,y2]=lin_fit2(x,y);
function [x1,y1,x2,y2]=lin_fit2(x,y);

%global s0
n=length(x(:));
mx=mean(x(:));
my=mean(y(:));
sxx=sum((x(:)-mx).^2)/(n-1);
syy=sum((y(:)-my).^2)/(n-1);
sxy=sum((x(:)-mx).*(y(:)-my))/(n-1);
d=1; % equal error variances, orthogonal regression
% best-fit slope:
s0=(syy-d*sxx+sqrt((syy-d*sxx)^2+4*d*sxy^2))/2/sxy;
% best-fit intercept:
s1=my-s0*mx;

xmin=min(x(:));
xmax=max(x(:));
ymin=min(y(:));
ymax=max(y(:));

y1=s1+xmin*s0;
y2=s1+xmax*s0;
x1=(y1-s1)/s0;
x2=(y2-s1)/s0;

%x1=(y1m-s1)/s0;
%x2=(y2m-s1)/s0;

if x1 < xmin 
 x1=xmin; 
end
if x2 > xmax
 x2=xmax;
end

if s0>0
 if y1 < ymin  
  x1=(ymin-s1)/s0;
  y1=ymin;
 end
 if y2 > ymax
  x2=(ymax-s1)/s0;
  y2=ymax;
 end
else
 if y1 > ymax 
  x1=(ymax-s1)/s0;
  y1=ymax;
 end
 if y2 < ymin
  x2=(ymin-s1)/s0;
  y2=ymin;
 end
end

return
