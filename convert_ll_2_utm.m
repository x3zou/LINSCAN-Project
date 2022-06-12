clear
clc
f = 'sc2019_hash_ABCD_so.focmec.scedc.txt';
% f = 'sc.focmec.scedc.txt';
fid = fopen(f);
C = textscan(fid,'%d %d %d %d %d %f %d %f %f %f %f %f %f %f %f %f %f %f %d %f %s', 'delimiter', '\n');
fclose(fid);

wid=length(C(1,:));
N=cell2mat(C(:,8:10));    % coordinates
FM=cell2mat(C(:,11:14)); % focal mech (M strike dip rake)
yr=cell2mat(C(:,1));     % year
mo=cell2mat(C(:,2));     % month
L=C{:,wid};              % quality flag

ind=find(yr< 2020);

% regional
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

x=N(ind,2);
y=N(ind,1);
z=-N(ind,3);
FM=FM(ind,1:4);

[xx,yy]=utm2ll(x,y,0,1);
xy=horzcat(xx-xo,yy-yo);
xy=xy*1e-3;
data=horzcat(yr,mo,N,FM);
data_converted=horzcat(yr,mo,(xx-xo)*1e-3,(yy-yo)*1e-3,z,FM);


