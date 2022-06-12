%load fm catalog and create sub-files
%Xiaoyu Zou, 6/2/2022
% ### column description of focmecall_SC_01_25_2021 file.
% 1 year
% 2 month
% 3 day
% 4 hour
% 5 minute
% 6 second
% 7 event ID
% 8 latitude
% 9 longitude
% 10 depth
% 11 magnitude
% 12 strike
% 13 dip
% 14 rake
% 15 1st nodal plane uncertainty
% 16 2nd nodal plane uncertainty
% 17 number of polarities
% 18 polarity misfit
% 19 number of S/P amplitude ratio
% 20 average log10(S/P) amplitude ratio misfit
% 21 focal mechanism quality
% 22 probability of solution close to real solution
% 23 azimuth gap
% 24 take-off angle gap
% 25 polarity misfit
% 26 station distribution ratio (STDR)
clear
clc
data=importdata("focmecall_SC_01_25_2021.txt");
data=data.textdata;
data(:,21)=[];
data=str2double(data);
data1=data(1:616483,:);%1981 to Jul 3 2019
data1=[data1(:,8) data1(:,9) data1(:,11) data1(:,12) data1(:,13) data1(:,14)];%xyMSDR
data2=data(616484:end,:);%Starting from Jul 4 2019
data2=[data2(:,8) data2(:,9) data2(:,11) data2(:,12) data2(:,13) data2(:,14)];

%local origin
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

% %Bounds for Ridgecrest Region: Easting[-60, 25], Northing[0,75]
[xx1,yy1]=utm2ll(data1(:,2),data1(:,1),0,1);%switched the position of latitude and longitude
[xx2,yy2]=utm2ll(data2(:,2),data2(:,1),0,1);
xy1=horzcat(xx1-xo,yy1-yo);
xy2=horzcat(xx2-xo,yy2-yo);
data1=[(xx1-xo)*1e-3 (yy1-yo)*1e-3 data1(:,3) data1(:,4) data1(:,5) data1(:,6)];%easting, northing, MSDR
data2=[(xx2-xo)*1e-3 (yy2-yo)*1e-3 data2(:,3) data2(:,4) data2(:,5) data2(:,6)];
d1=zeros(1,6);
d2=zeros(1,6);
for i=1:length(data1)
    if (data1(i,1)>-60) && (data1(i,1)<25) && (data1(i,2)>0) && (data1(i,2)<75)
        d1=vertcat(d1,data1(i,:));
    end
end
d1(1,:)=[];
for i=1:length(data2)
    if (data2(i,1)>-60) && (data2(i,1)<25) && (data2(i,2)>0) && (data2(i,2)<75)
        d2=vertcat(d2,data2(i,:));
    end
end
d2(1,:)=[];

d1xy=[d1(:,1),d1(:,2)];
d2xy=[d2(:,1),d2(:,2)];
writematrix(d1xy,'date1input.txt')
writematrix(d2xy,'date2input.txt')
writematrix(d1,'date1MSDR.txt')
writematrix(d2,'date2MSDR.txt')




