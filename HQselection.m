%This is the algorithm for selecting high-quality clusters
%Xiaoyu Zou, x3zou@ucsd.edu,07/07/2022
clear
clc
list=importdata("d2listeps=sqrt8_2.txt");%lin_scan output list,change with iteration
%real=importdata("eps=sqrt(16)_input.txt");
real=importdata("d2eps=sqrt8_input_2.txt");%lin_scan input,change with iteration
points=importdata("d2eps=sqrt8_2.txt");%lin_scan output,change with iteration
catalog=importdata("date2MSDR.txt");%complete earthquake catalog
input=importdata("date2input.txt");%earthquake catalog xy data
points_selected=zeros(1,2);
points_experiment=zeros(1,2);
points_added=zeros(1,2);
clusters_saved=zeros(1,3);
kldwhat=zeros(1,1);
dipP1test=zeros(1,1);
dipP2test=zeros(1,1);
dipT1test=zeros(1,1);
dipT2test=zeros(1,1);
z=56;%!!!please rememeber you may need to change this value everytime you change the importdata!!!
ends=zeros(1,4);
for i=min(list):max(list)
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
        %x=rmoutliers(x,'percentiles',[10 90]);
        x=rmoutliers(x,'quartiles');
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        pts=[x(:,1) x(:,2) 0*x(:,1)];
        L=sqrt((x2-x1)^2+(y2-y1)^2);
        [d] = point_to_line(pts, v1, v2);
        mean_dist = mean(abs(d))/L;
    else
        continue
    end
    if isnan(L)
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
% m=15;
% n=15;
% ha = tight_subplot(m,n,[0.03 0.03],[0.1 0.01],[0.01 0.01]);
k=0;
for i=min(list):max(list)
    pointsdr=zeros(1,6);
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
%         x=rmoutliers(x,'percentiles',[10 90]);
        x=rmoutliers(x,'quartiles');
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        L=sqrt((x1-x2)^2+(y1-y2)^2);
        slope=atan((y1-y2)/(x1-x2));
    else
        continue
    end
    if isnan(slope)
        continue
    end
    if L>20
        continue
    end
    %eliminate the marginal effect (tentatively)
    if (abs(slope)<5*pi/180 || abs(slope)>pi/2.1) && (mean(x(:,1))> (max(input(:,1)-2)) || mean(x(:,1))< (min(input(:,1)+2)) || mean(x(:,2))> (max(input(:,2))-2) || mean(x(:,2))<(min(input(:,2))+2))
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
    if kld < 1 && (abs(cc) > 0.7 || mean_dist_2 <0.08)%for uniform distribution
        for c = 1:size(pp,1)
            g=pp(c,1);
            h=pp(c,2);
            ind=find(input(:,1)==g & input(:,2)==h);
            if size(ind,1) > 1
                g=repmat(g,size(ind,1),1);
                h=repmat(h,size(ind,1),1);
            end
            M=catalog(ind,3);
            S=catalog(ind,4);
            D=catalog(ind,5);
            R=catalog(ind,6);
            xySDR=horzcat(g,h,M,S,D,R);
            pointsdr=vertcat(pointsdr,xySDR);
        end
        pointsdr(1,:)=[];
        fm=pointsdr(:,3:6);
        arr=[pointsdr(:,1) pointsdr(:,2) fm];
        deg2rad=pi/180.;
        comp=zeros(1,6);
        sm=[];
        num_events=size(arr(:,1),1);
        for k=1:num_events
            [fm]=sdr2mij(arr(k,4:6),deg2rad);
            Mo=10^(1.5*(arr(k,3)+10.7));
            fm=fm*Mo;
            comp=comp+fm;
            % this is to estimate variability:
            %[s1,d1,s2,d2] = mij2sd(fm);
            %sm=[sm min([s1 s2])];
        end
        %stdm=std(sm);
        Mo=sqrt(sum(comp(:).^2)/2);
        comp=comp/Mo;
        [S1,D1,R1,S2,D2,R2]=mij2sdr2(comp);
        [azP1, dipP1, azT1, dipT1] = PTaxis(S1,D1,R1);
        [azP2, dipP2, azT2, dipT2] = PTaxis(S2,D2,R2);
        dipT1test=vertcat(dipT1test,dipT1);
        dipT2test=vertcat(dipT2test,dipT2);
        dipP1test=vertcat(dipT1test,dipP1);
        dipP2test=vertcat(dipT1test,dipP2);
        if dipP1<45*pi/180 && dipT1<45*pi/180 && dipP2 < 45*pi/180 && dipT2 < 45*pi/180 %for strike-slip fault
            k=k+1;
            z=z+1;
            j=k;
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
title("eps=sqrt(4_2)")
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

