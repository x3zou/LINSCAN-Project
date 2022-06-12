clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; % area bounds
Fs=14; % font size

deg2rad=pi/180.;

% local origin
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

% data file
f = 'sc.focmec.scedc.txt';
fid = fopen(f);
C = textscan(fid,'%d %d %d %d %d %f %d %f %f %f %f %f %f %f %f %f %f %f %d %f %s', 'delimiter', '\n');
fclose(fid);

wid=length(C(1,:));
N=cell2mat(C(:,8:10));    % coordinates
FM=cell2mat(C(:,12:14)); % focal mech (strike dip rake)
yr=cell2mat(C(:,1));     % year
mo=cell2mat(C(:,2));     % month
L=C{:,wid};              % quality flag

%ind=find((yr< 2019 |(yr==2019 & mo<7)) );
ind=find(yr< 2020);
%ind=find(strtrim(char(L)) == 'A');
%ind=find(strtrim(char(L)) == 'A' | strtrim(char(L)) == 'B');  % select high-quality data

x=N(ind,2);  % longtitude
y=N(ind,1);  % latitude
z=-N(ind,3); % depth (negative down)

FM=FM(ind,1:3); % fault attitude: strike, dip, rake

qn=0*x;
Q=strtrim(char(L(ind)));     
% assign numerical score to quality flag 
ind=find(Q == 'A');
qn(ind)=4;
ind=find(Q == 'B');
qn(ind)=3;
ind=find(Q == 'C');
qn(ind)=2;
ind=find(Q == 'D');
qn(ind)=1;

%return


% convert lat/lon to local Cartesian coordinates
[xx,yy]=utm2ll(x(:),y(:),0,1);

xu=reshape((xx-xo)*1e-3,size(x));
yu=reshape((yy-yo)*1e-3,size(x));

% line(xu,yu,'lineStyle','none','Marker','.','MarkerSize',1,'MarkerFaceColor',grey,'MarkerEdgeColor',grey), hold on

% select data from a given region
ind=find(xu < maxlon & xu > minlon & yu < maxlat & yu > minlat );

xu=xu(ind);
yu=yu(ind);
z=z(ind);
FM=FM(ind,1:3);

%plot3(xu,yu,z,'k.');

figure(1)
clf
%for i=1:length(xu)
for i=1:10 
 x1=xu(i); y1=yu(i); z1=z(i); 
 plot3(x1,y1,z1,'k.');
 xb= [x1-1 x1-1 x1+1 x1+1];
 yb= [y1-1 y1-1 y1+1 y1+1];
 zb= [z1-1 z1+1 z1+1 z1-1];

 fill3(xb,yb,zb,1,'EdgeColor','k','LineWidth',1,'FaceAlpha',0.5,'EdgeAlpha',1); hold on
end

axis equal;
