clear
list=importdata("02list_eps=sqrt2.txt");
%real=importdata("eps=sqrt(16)_input.txt");
real=importdata("realdata.txt");
points=importdata("correlation0.2iter_eps=sqrt2.txt");
points_selected=zeros(1,2);
points_experiment=zeros(1,2);
points_added=zeros(1,2);
clusters_saved=zeros(1,3);
kldwhat=zeros(1,1);
z=0;%!!!please rememeber you may need to change this value everytime you change the importdata!!!
ends=zeros(1,4);
for i=min(list):max(list)
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        pts=[x(:,1) x(:,2) 0*x(:,1)];
        L=sqrt((x2-x1)^2+(y2-y1)^2);
        [d] = point_to_line(pts, v1, v2);
        mean_dist = mean(abs(d))/L;
    else
        continue
    end
    if L>20
        continue
    end
    points_experiment=vertcat(points_experiment,x);
    %         if mean_dist > 0.01 %original value:0.05
    points_selected=vertcat(points_selected,x);
    %         end
end
points_selected(1,:)=[];
idx=~ismember(real,points_selected,'rows');
res=real(idx,:);
%writematrix(res,'selected_eps=sqrt10.txt')
%writematrix(points_selected,'points_selected_eps=sqrt14.txt')
m=15;
n=15;
ha = tight_subplot(m,n,[0.03 0.03],[0.1 0.01],[0.01 0.01]);
k=0;
for i=min(list):max(list)
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        L=sqrt((x1-x2)^2+(y1-y2)^2);
    else
        continue
    end
    if L>20
        continue
    end
    pts=[x(:,1) x(:,2) 0*x(:,1)];
    pts_1=[res(:,1) res(:,2) 0*res(:,1)];
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    [d] = point_to_line(pts, v1, v2);
    %[d1]= point_to_line(pts_1,v1,v2);
    mean_dist=mean(abs(d))/L;
    %norm_dist = abs(d1)./L;
    %if mean_dist > 0.01 %original value: 0.05
    [ind,~]=point_to_line2(pts_1,v1,v2,0.1); %original value:0.05
    add=[res(ind,1),res(ind,2)];
    pp=vertcat(add,x);
    cc=pearson(pp(:,1),pp(:,2));
    pts_2=[pp(:,1) pp(:,2) 0*pp(:,1)];
    [d2]=point_to_line(pts_2,v1,v2);
    mean_dist_2=mean(abs(d2))/L;
    kld=getkld(pp);
    kldwhat=vertcat(kldwhat,kld);
    if kld < 1
        if abs(cc) > 0.7 || mean_dist_2 <0.07 %original value: 0.8, 0.05
            kld
            k=k+1;
            z=z+1;
            j=k;
            axes(ha(j));
            %                 line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
            %                 line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b')
            %                 line([x1 x2],[y1 y2],'lineStyle','-','Color','b','LineWidth',1)
            ends=vertcat(ends,[x1 y1 x2 y2]);
            %                 axis equal
            %                 axis manual
            %                 line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
            %                 line(points(:,1),points(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
            %                 line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
            %                 line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','b','MarkerEdgeColor','b')
            num=size(pp,1);
            length=sqrt((x1-x2)^2+(y1-y2)^2);
            %                 text((x1+x2)/2,(y1+y2)/2,sprintf('n=%i,d=%f,r=%f,l=%f',num,mean_dist_2,cc,length),'Color','r','FontWeight','bold');
            points_added = vertcat(points_added,pp);
            clusters_added=[z*ones(size(pp,1),1) pp(:,1) pp(:,2)];
            clusters_saved=vertcat(clusters_saved,clusters_added);
            %points_added(1,:)=[];
            %points_total=vertcat(pp,x);
            %clusters_added=[z*ones(size(points_added,1),1) points_added(:,1) points_added(:,2)];
            %clusters_saved=vertcat(clusters_saved,clusters_added);
        end
    end
end
points_added(1,:)=[];
%clusters_added(1,:)=[];
clusters_saved(1,:)=[];
ends(1,:)=[];
pttest=[add(:,1) add(:,2) 0*add(:,1)];
[dtest]=point_to_line(pttest,v1,v2);
test_normdist=abs(dtest)./L;

figure()
line(points(:,1),points(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(points_added(:,1),points_added(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
legend("LINSCAN PICKUP","HIGH QUALITY")
title("eps=sqrt(16)")
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

