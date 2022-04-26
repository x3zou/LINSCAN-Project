%clear all;
close all;

%load cluster6.txt;
%load 0.2iter_eps=sqrt5.mat
x=c75;
data=importdata('realdata.txt');
clf

scatter(x(:,1),x(:,2),'k'), hold on
[x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
line([x1 x2],[y1 y2]), hold on
v1=[x1 y1 0];  v2=[x2 y2 0];
pts=[x(:,1) x(:,2) 0*x(:,1)];
L=sqrt((x2-x1)^2+(y2-y1)^2);
[d] = point_to_line(pts, v1, v2);

mean_dist = mean(abs(d))/L;
text((x1+x2)/2,(y1+y2)/2,sprintf('mean = %6.3f \n',mean_dist-0.05));

axis manual;
scatter(data(:,1),data(:,2),'r')
scatter(x(:,1),x(:,2),'k')
line([x1 x2],[y1 y2])
title("cluster1")


%saveas(gcf,'cluster89.png')
return

function [d] = point_to_line(pt, v1, v2)
% function [d] = point_to_line(pt, v1, v2)
% returns distances from points to a line 

 L = norm(v1 - v2);

 v1 = repmat(v1,size(pt,1),1); % 1st end of line point
 v2 = repmat(v2,size(pt,1),1); % 2nd end of line point
 a = v1 - v2;
 b = pt - v2;
 c = pt - v1;
 d = sqrt(sum(cross(a,b,2).^2,2))/L; %distance to line

end


