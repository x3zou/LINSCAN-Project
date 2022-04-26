clear
real=importdata('realdata.txt');
xtest=importdata("kldtesting.txt");
% x1=importdata("sqrt2_1_41.txt");
% x2=importdata("sqrt4_42_173.txt");
% x3=importdata("sqrt8_174_231.txt");
% x4=importdata("sqrt16_232_235.txt");
% endtest=importdata("end_test.txt");
% end1=importdata("sqrt2_end.txt");
% end2=importdata("sqrt4_end.txt");
% end3=importdata("sqrt8_end.txt");
% end4=importdata("sqrt16_end.txt");
% endpoints=vertcat(end1,end2,end3,end4);
endpoints=importdata("kldends.txt");
% all=vertcat(x1,x2,x3,x4);
all=xtest;
new_endpoints=zeros(1,4);
m=6;
n=6;
j=0;
ha = tight_subplot(m,n,[0.03 0.03],[0.1 0.01],[0.01 0.01]);
for i=1:m*n
    axes(ha(i));
    j=j+1;
%     if j==3 || j==5 || j==8 || j==17 ||j==48 || j==51 || j==54 || j==61|| j==68 ||j==71|| j==72
%         j=j+1;
%     end
    pp=horzcat(all(find(all(:,1)==j),2),all(find(all(:,1)==j),3));
    kld=getkld(pp);
    x1=endpoints(j,1);
    x2=endpoints(j,3);
    y1=endpoints(j,2);
    y2=endpoints(j,4);
    v1=[x1 y1 0];  v2=[x2 y2 0];
    pts=[pp(:,1) pp(:,2) 0*pp(:,1)];
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    [d] = point_to_line(pts, v1, v2);
    mean_dist=mean(abs(d))/L;
    line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
    line([x1 x2],[y1 y2],'lineStyle','-','Color','b','LineWidth',1)
    axis equal
    axis manual
    line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
    line(all(:,2),all(:,3),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
    line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k') 
    num=size(pp,1);
    cc=pearson(pp(:,1),pp(:,2));
    no=i;
    text((x1+x2)/2,(y1+y2)/2,sprintf('n=%i,d=%f,r=%f,l=%f,kld=%f',num,mean_dist,cc,L,kld),'Color','r','FontWeight','bold','Fontsize',8);
    newendpoints=[x1 y1 x2 y2];
    new_endpoints=vertcat(new_endpoints,newendpoints);
end
new_endpoints(1,:)=[];
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

ind=find((d1+d2)/L<D+sqrt(1+D^2)); % index to nearby points

end

