%Plot the distribution of fault traces and focal mechanisms. 
%Also calculate the input for dihedral angle calculation
%Xiaoyu Zou, x3zou@ucsd.edu, 7/31/2022
clear
clc
minlon = 0; maxlon = 60; minlat = -200; maxlat = 2; %scope limit for fault traces plotting. Unit: UTM
%list=importdata("02list_eps=sqrt4_test.txt");
ends=importdata("ttnewends.txt");
%points=importdata("correlation0.2test.txt");%linscan-scanned data
points2=importdata("ttnewclusters.txt");%high quality clusters
input=importdata("ttinput.txt");%earthquake catalog xy data
catalog=importdata("ttMSDR.txt");%complete earthquake catalog
totalcomp=importdata('ttnewcomp.txt');%composite focal mechanism for each cluster
pointsdr=zeros(1,7);
newpointsdr=zeros(1,6);
newends=zeros(1,4);
ERROR=zeros(1,1);
SDRSDR=zeros(1,6);
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
    M=catalog(ind,4);
    S=catalog(ind,5);
    D=catalog(ind,6);
    R=catalog(ind,7);
    xySDR=horzcat(n,x,y,M,S,D,R);
    pointsdr=vertcat(pointsdr,xySDR);
end
pointsdr(1,:)=[];
j=0;
% line(input(:,1),input(:,2),'lineStyle','none','Marker','.','MarkerSize',6,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
% hold on
for i=1:max(points2(:,1))
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
%         [s1,d1,s2,d2] = mij2sd(fm);
%         sm=[sm min([s1 s2])];
    end
%     stdm=std(sm);
    P=zeros(1,2);
    P(1)=(y1-y2)/(x1-x2);
    P(2)=y1-P(1)*x1;
    err=err_2D_data(j,pp(:,1),pp(:,2),P);
    Mo=sqrt(sum(comp(:).^2)/2);
    comp=comp/Mo;
    comp=totalcomp(j,:);
    [S1,D1,R1,S2,D2,R2]=mij2sdr2(comp);
    [azP1, dipP1, azT1, dipT1] = PTaxis(S1,D1,R1);
    [azP2, dipP2, azT2, dipT2] = PTaxis(S2,D2,R2);
    if dipP1<=40*pi/180 && dipT1<=40*pi/180 && dipP2 <= 40*pi/180 && dipT2 <= 40*pi/180 %for strike-slip fault
     SDRSDR=vertcat(SDRSDR,[S1 D1 R1 S2 D2 R2]);
     ERROR=vertcat(ERROR,err);
     newends=vertcat(newends,[x1,y1,x2,y2]);
    else
        continue
    end
%     v1=[x1 y1 0];  v2=[x2 y2 0];
%     pts=[pp(:,1) pp(:,2) 0*pp(:,1)];
%     L=sqrt((x2-x1)^2+(y2-y1)^2);
%     [d] = point_to_line(pts, v1, v2);
%     mean_dist=mean(abs(d))/L;
%     num=size(pp,1);
%     cc=pearson(pp(:,1),pp(:,2));
%     no=i;
%     M = [comp(1) comp(4) comp(5); comp(4) comp(2) comp(6); comp(5) comp(6) comp(3)];
%     %this is to determine where to plot the beachball
%     slope=(y1-y2)/(x1-x2);
%     if slope>0
%         xm=(x1+x2)/2+L/3.5; ym=(y1+y2)/2-L/5;
%     else
%         xm=(x1+x2)/2-L/3.5; ym=(y1+y2)/2-L/5;
%     end
%     line(pp(:,1),pp(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k')
%     focalmech(comp, xm, ym, 3.5, 'b','dc')%plot focal mechanism distribution
%     hold on
%     axis equal
end
% legend("background","clusters")
ERROR(1,:)=[];
SDRSDR(1,:)=[];
newends(1,:)=[];

% input for dihedral angle calculation
ttdihedral_input=horzcat(newends,SDRSDR,ERROR);
writematrix(ttdihedral_input,"ttdihedral_newinput.txt")

%plot california faults
orig=[-117.5 35.5]; % <— change to your local origin
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
fault_tr = 1; % flag for plotting the fault trace
col=[1 0 1]; % line color
if fault_tr==1
    load ca_faults_ll.dat;
    faults=ca_faults_ll;
    load ca_dim_ll.dat;
    dim=ca_dim_ll;
    [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
    faults(:,1)=(faults(:,1)-xo)*1e-3;
    faults(:,2)=(faults(:,2)-yo)*1e-3;
    [trac]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],col);
end
xlabel("Eastings(km)")
ylabel("Northings(km)")
title("trans-pressional ")
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