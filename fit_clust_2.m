clear all;
close all;

load cluster6.txt;
x=cluster6;
clf

plot(x(:,1),x(:,2),'bo'), hold on
[x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
plot([x1 x2],[y1 y2],'m-'), hold on
v1=[x1 y1 0];  v2=[x2 y2 0];
pts=[x(:,1) x(:,2) 0*x(:,1)];
L=sqrt((x2-x1)^2+(y2-y1)^2);
%[d] = point_to_line(pts, v1, v2);
%mean_dist = mean(abs(d))/L;
%fprintf('mean = %6.3f \n',mean_dist);

mean_dist=0.1;
[ind,d] = point_to_line2(pts, v1, v2, mean_dist);

plot(x(ind,1),x(ind,2),'k+'), hold on

axis equal;


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


function [ind,d] = point_to_line2(pt, v1, v2, D)
% function [ind,d] = point_to_line2(pt, v1, v2, D)
% return points within relative distance D to a finite line segment

 L = norm(v1 - v2);

 v1 = repmat(v1,size(pt,1),1); % 1st end of line point
 v2 = repmat(v2,size(pt,1),1); % 2nd end of line point
 a = v1 - v2;
 b = pt - v2;
 c = pt - v1;
 d1= sqrt(sum(c.*c,2)); % square of distance to one end of the line segment
 d2= sqrt(sum(b.*b,2)); % square of distance to the other end of the line segment
 d = sqrt(sum(cross(a,b,2).^2,2))/L; %distance to line

ind=find(d/L < D & (d1+d2)/L < 1+3*D); % index to nearby points

end

