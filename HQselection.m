%This is the algorithm for selecting high-quality clusters
%Xiaoyu Zou, x3zou@ucsd.edu,07/31/2022
clear
clc
list=importdata("ttlisteps=sqrt16.txt");%lin_scan output list,change with iteration
%real=importdata("eps=sqrt(16)_input.txt");
real=importdata("tteps=16_input.txt");%lin_scan input,change with iteration
points=importdata("tteps=sqrt16.txt");%lin_scan output,change with iteration
catalog=importdata("ttMSDR.txt");%complete earthquake catalog
input=importdata("ttinput.txt");%earthquake catalog xy data
points_selected=zeros(1,2);
points_experiment=zeros(1,2);
points_added=zeros(1,2);
points_added2=zeros(1,2);
clusters_saved=zeros(1,3);
clusters_saved2=zeros(1,3);
clusters_saved3=zeros(1,3);
clusters_saved4=zeros(1,3);
totalcomp=zeros(1,6);
kldwhat=zeros(1,1);
dipP1test=zeros(1,1);
dipP2test=zeros(1,1);
dipT1test=zeros(1,1);
dipT2test=zeros(1,1);
testingdip=zeros(1,4);
z=0;%don't change this value
ends=zeros(1,4);
ends2=zeros(1,4);
ends3=zeros(1,4);
for i=min(list):max(list)
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
        x=rmoutliers(x,'quartiles');
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        pts=[x(:,1) x(:,2) 0*x(:,1)];
        L=sqrt((x2-x1)^2+(y2-y1)^2);
        [d] = point_to_line(pts, v1, v2);
        mean_dist = mean(abs(d))/L;%mean normalized distance
    else
        continue
    end
    if isempty(x)
        continue
    end
    if size(x,1)<10
        continue
    end
    if isnan(L)
        continue
    end
    if L>15
        continue
    end
    points_experiment=vertcat(points_experiment,x);
    points_selected=vertcat(points_selected,x);
end
points_selected(1,:)=[];
idx=~ismember(real,points_selected,'rows');
res=real(idx,:);
for i=min(list):max(list)
    pointsdr=zeros(1,6);
    a=double(isempty(points(find(list==i),:)));
    if a==0
        x=points(find(list==i),:);
        x=rmoutliers(x,'quartiles');
        [x1,y1,x2,y2]=lin_fit2(x(:,1),x(:,2));
        v1=[x1 y1 0];  v2=[x2 y2 0];
        [d] = point_to_line(pts, v1, v2);
        L=sqrt((x1-x2)^2+(y1-y2)^2);
        slope=atan((y1-y2)/(x1-x2));
    else
        continue
    end
    if isempty(x)
        continue
    end
    if size(x,1)<10
        continue
    end
    if isnan(slope)
        continue
    end
    if L>15
        continue
    end
    %eliminate the marginal effect (tentatively)
    if (abs(slope)<5*pi/180 || abs(slope)>pi/2.1) && (mean(x(:,1))> (max(input(:,1)-2)) || mean(x(:,1))< (min(input(:,1)+2)) || mean(x(:,2))> (max(input(:,2))-2) || mean(x(:,2))<(min(input(:,2))+2))
        continue
    end
    pts=[x(:,1) x(:,2) 0*x(:,1)];
    pts_1=[res(:,1) res(:,2) 0*res(:,1)];
    L=sqrt((x2-x1)^2+(y2-y1)^2);
    mean_dist=mean(abs(d))/L;
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
            M=catalog(ind,4);
            S=catalog(ind,5);
            D=catalog(ind,6);
            R=catalog(ind,7);
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
        if dipP1<=40*pi/180 && dipT1<=40*pi/180 && dipP2 <= 40*pi/180 && dipT2 <= 40*pi/180 %for strike-slip fault
            k=k+1;
            z=z+1;
            j=k;
            ends=vertcat(ends,[x1 y1 x2 y2]);
            num=size(pp,1);
            length=sqrt((x1-x2)^2+(y1-y2)^2);
            points_added = vertcat(points_added,pp);
            clusters_added=[z*ones(size(pp,1),1) pp(:,1) pp(:,2)];
            x_2=[clusters_added(:,2) clusters_added(:,3)];
            idx=~ismember(real,x_2,'rows');
            res2=real(idx,:);
            pts_3=[res2(:,1) res2(:,2) 0*res2(:,1)];
            [ind,~]=point_to_line2(pts_3,v1,v2,0.1); %original value:0.05
            add=[res2(ind,1),res2(ind,2)];
            pp2=vertcat(add,x_2);
            clusters_added2=[z*ones(size(pp2,1),1) pp2(:,1) pp2(:,2)];
            clusters_saved=vertcat(clusters_saved,clusters_added);
            clusters_saved2=vertcat(clusters_saved2,clusters_added2);
        end
    end
end
clusters_saved2(1,:)=[];
clusters_saved(1,:)=[];
z=96;%change z here
ends(1,:)=[];
for i=1:size(ends,1)
    pp=[clusters_saved2(find(clusters_saved2(:,1)==i),2),clusters_saved2(find(clusters_saved2(:,1)==i),3)];
    x1=ends(i,1);
    y1=ends(i,2);
    x2=ends(i,3);
    y2=ends(i,4);
    v1=[x1 y1 0];
    v2=[x2,y2 0];
    length=sqrt((x1-x2)^2+(y1-y2)^2);
    pts_3=[pp(:,1) pp(:,2) 0*pp(:,1)];
    [d3]=point_to_line(pts_3,v1,v2);
    pp=pp(find(abs(d3)/length<0.1),:);%remove points too far from best fitting
    [x1,y1,x2,y2]=lin_fit2(pp(:,1),pp(:,2));%fit a new line with modified cluster
    cc=pearson(pp(:,1),pp(:,2));
    mean_dist_3=mean(abs(d3))/L;
    kld=getkld(pp);
    if kld <= 1 && (abs(cc) >= 0.7)%for uniform distribution
        %     if kld<1
        for c = 1:size(pp,1)
            g=pp(c,1);
            h=pp(c,2);
            ind=find(input(:,1)==g & input(:,2)==h);
            if size(ind,1) > 1
                g=repmat(g,size(ind,1),1);
                h=repmat(h,size(ind,1),1);
            end
            M=catalog(ind,4);
            S=catalog(ind,5);
            D=catalog(ind,6);
            R=catalog(ind,7);
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
        [azP1, dipP1, azT1, dipT1] = PTaxis(S1,D1,R1);
        [azP2, dipP2, azT2, dipT2] = PTaxis(S2,D2,R2);
        if (dipP1<=40*pi/180 && dipT1<=40*pi/180 && dipP2 <= 40*pi/180 && dipT2 <= 40*pi/180) && (abs(S1-strike)<=20 ||abs(S2-strike)<=20 || abs(S1-180-strike)<=20 || abs(S2-180-strike)<=20)
                z=z+1;
                ends2=vertcat(ends2,[x1 y1 x2 y2]);
                num=size(pp,1);
%                 dipinfo=[dipP1*180/pi dipP2*180/pi dipT1*180/pi dipT2*180/pi];
%                 testingdip=[testingdip;dipinfo];
                totalcomp=[totalcomp;comp];
                length=sqrt((x1-x2)^2+(y1-y2)^2);
                points_added2 = vertcat(points_added2,pp);
                clusters_added=[z*ones(size(pp,1),1) pp(:,1) pp(:,2)];
                clusters_saved3=vertcat(clusters_saved3,clusters_added);
        end
    end
end
points_added(1,:)=[];
points_added2(1,:)=[];
clusters_saved3(1,:)=[];
number=0;
testingdip(1,:)=[];
totalcomp(1,:)=[];
% k=65;%change the value here!
% for i=1:z
%     clusters_saved3(find(clusters_saved3(:,1)==i),1)=k+i;
% end
ends2(1,:)=[];

% %for selecting strike-slip clusters:(don't put this output into iteration
% %steps, use clusters_saved3)
% z=0;%remember to change it!!
% for i=1:size(ends2,1)
%     pp=[clusters_saved3(find(clusters_saved3(:,1)==i),2),clusters_saved3(find(clusters_saved3(:,1)==i),3)];
%     x1=ends2(i,1);
%     y1=ends2(i,2);
%     x2=ends2(i,3);
%     y2=ends2(i,4);
%         for c = 1:size(pp,1)
%             g=pp(c,1);
%             h=pp(c,2);
%             ind=find(input(:,1)==g & input(:,2)==h);
%             if size(ind,1) > 1
%                 g=repmat(g,size(ind,1),1);
%                 h=repmat(h,size(ind,1),1);
%             end
%             M=catalog(ind,4);
%             S=catalog(ind,5);
%             D=catalog(ind,6);
%             R=catalog(ind,7);
%             xySDR=horzcat(g,h,M,S,D,R);
%             pointsdr=vertcat(pointsdr,xySDR);
%         end
%         pointsdr(1,:)=[];
%         fm=pointsdr(:,3:6);
%         arr=[pointsdr(:,1) pointsdr(:,2) fm];
%         deg2rad=pi/180.;
%         comp=zeros(1,6);
%         sm=[];
%         num_events=size(arr(:,1),1);
%         for k=1:num_events
%             [fm]=sdr2mij(arr(k,4:6),deg2rad);
%             Mo=10^(1.5*(arr(k,3)+10.7));
%             fm=fm*Mo;
%             comp=comp+fm;
%             % this is to estimate variability:
%             %[s1,d1,s2,d2] = mij2sd(fm);
%             %sm=[sm min([s1 s2])];
%         end
%         %stdm=std(sm);
%         Mo=sqrt(sum(comp(:).^2)/2);
%         comp=comp/Mo;
%         [S1,D1,R1,S2,D2,R2]=mij2sdr2(comp);
%         [azP1, dipP1, azT1, dipT1] = PTaxis(S1,D1,R1);
%         [azP2, dipP2, azT2, dipT2] = PTaxis(S2,D2,R2);
%         if dipP1<44*pi/180 && dipT1<44*pi/180 && dipP2 < 44*pi/180 && dipT2 < 44*pi/180 %for strike-slip fault
% %             if abs(S1-strike)<=30 ||abs(S2-strike)<=30 || abs(S1-180-strike)<=30 || abs(S2-180-strike)<=30
%                 z=z+1;
%                 ends3=vertcat(ends3,[x1 y1 x2 y2]);
%                 num=size(pp,1);
%                 length=sqrt((x1-x2)^2+(y1-y2)^2);
%                 clusters_added=[z*ones(size(pp,1),1) pp(:,1) pp(:,2)];
%                 clusters_saved4=vertcat(clusters_saved4,clusters_added);
% %             end
%         end
% end
% ends3(1,:)=[];
% clusters_saved4(1,:)=[];
figure()
line(points(:,1),points(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(clusters_saved3(:,2),clusters_saved3(:,3),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
legend("LINSCAN PICKUP","HIGH QUALITY")
title("eps=sqrt(8)")
% 
writematrix(ends2,'testingends.txt')
writematrix(clusters_saved3,'testingclusters.txt')
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

