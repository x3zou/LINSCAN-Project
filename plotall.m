%This is the algorithms for plotting all indivicual clusters and focal
%mechanisms
%Xiaoyu Zou, x3zou@ucsd.edu, 7/11/2022
clear
clc
%list=importdata("02list_eps=sqrt4_test.txt");
ends=importdata("d2ends.txt");
%points=importdata("correlation0.2test.txt");%linscan-scanned data
points2=importdata("d2clusters.txt");%high quality clusters
input=importdata("date2input.txt");%earthquake catalog xy data
catalog=importdata("date2MSDR.txt");%complete earthquake catalog
pointsdr=zeros(1,7);
newpointsdr=zeros(1,6);
%Unify the sequence of endpoints
y1=ends(:,2);
y2=ends(:,4);
x1=ends(:,1);
x2=ends(:,3);
for i=1:length(ends)
    if y1(i)>y2(i)
        y1(i)=ends(i,4);
        y2(i)=ends(i,2);
        x1(i)=ends(i,3);
        x2(i)=ends(i,1);
    end
end
ends(:,2)=y1;
ends(:,4)=y2;
ends(:,1)=x1;
ends(:,3)=x2;

%Organize strike-dip-rake information
for i = 1:size(points2,1)
    n=points2(i,1);
    x=points2(i,2);
    y=points2(i,3);
    ind=find(input(:,1)==x & input(:,2)==y);
    if size(ind,1) > 1
        x=repmat(x,size(ind,1),1);
        y=repmat(y,size(ind,1),1);
        n=repmat(n,size(ind,1),1);
    end
    M=catalog(ind,3);
    S=catalog(ind,4);
    D=catalog(ind,5);
    R=catalog(ind,6);
    xySDR=horzcat(n,x,y,M,S,D,R);
    pointsdr=vertcat(pointsdr,xySDR);
end
pointsdr(1,:)=[];

%plotall
m=3;
n=3;
for p=1:11
    j=0+(p-1)*m*n;
    figure()
    ha = tight_subplot(m,n,[0.03 0.03],[0.1 0.01],[0.01 0.01]);
    for i=1:m*n
        axes(ha(i));
        j=j+1;
        pp=horzcat(pointsdr(find(pointsdr(:,1)==j),2),pointsdr(find(pointsdr(:,1)==j),3),pointsdr(find(pointsdr(:,1)==j),4),pointsdr(find(pointsdr(:,1)==j),5),pointsdr(find(pointsdr(:,1)==j),6),pointsdr(find(pointsdr(:,1)==j),7));
        x1=ends(j,1);
        x2=ends(j,3);
        y1=ends(j,2);
        y2=ends(j,4);
        fm=pp(:,3:6);
        arr=[pp(:,1) pp(:,2) fm];
        deg2rad=pi/180.;
        comp=zeros(1,6);
        sm=[];
        num_events=length(arr(:,1));
        for k=1:num_events
            [fm]=sdr2mij(arr(k,4:6),deg2rad);
            Mo=10^(1.5*(arr(k,3)+10.7));
            fm=fm*Mo;
            comp=comp+fm;
            % this is to estimate variability:
            [s1,d1,s2,d2] = mij2sd(fm);
            sm=[sm min([s1 s2])];
        end
        stdm=std(sm);
        Mo=sqrt(sum(comp(:).^2)/2);
        comp=comp/Mo;

        %estimate the composite strike
        strike=atan2(y2-y1,x2-x1)*180/pi;
        if strike<0
            strike=strike+180;
        end
        if strike <90
            strike = 90-strike;
        end
        if strike >90
            strike = 450-strike;
        end
        [S1,D1,R1,S2,D2,R2]=mij2sdr2(comp);
        [azP1, dipP1, azT1, dipT1] = PTaxis(S1,D1,R1);
        [azP2, dipP2, azT2, dipT2] = PTaxis(S2,D2,R2);
        if abs(strike-S1) > abs(strike-S2)
            dipP=dipP2;
        else
            dipP=dipP1;
        end

        %plot clusters and lines
        v1=[x1 y1 0];  v2=[x2 y2 0];
        pts=[pp(:,1) pp(:,2) 0*pp(:,1)];
        L=sqrt((x2-x1)^2+(y2-y1)^2);
        [d] = point_to_line(pts, v1, v2);
        mean_dist=mean(abs(d))/L;
        line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
        line([x1 x2],[y1 y2],'lineStyle','-','Color','b','LineWidth',1)
        axis equal
        axis manual
        line(input(:,1),input(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
        line(pointsdr(:,2),pointsdr(:,3),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0])
        line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
        num=size(pp,1);
        cc=pearson(pp(:,1),pp(:,2));
        no=i;

        M = [comp(1) comp(4) comp(5); comp(4) comp(2) comp(6); comp(5) comp(6) comp(3)];
        %this is to determine where to plot the beachball
        slope=(y1-y2)/(x1-x2);
        if slope>0
            xm=(x1+x2)/2+L/3.5; ym=(y1+y2)/2-L/5;
        else
            xm=(x1+x2)/2-L/3.5; ym=(y1+y2)/2-L/5;
        end
        focalmech(comp, xm, ym, L/7, 'b','dc')
        text((x1+x2)/2,(y1+y2)/2,sprintf('n=%i,d=%f,cc=%f,dipP=%f',num,mean_dist,cc,dipP),'Color','r','FontWeight','bold','Fontsize',8);
    end
end
newpointsdr=vertcat(newpointsdr,pp);
newpointsdr(1,:)=[];
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



