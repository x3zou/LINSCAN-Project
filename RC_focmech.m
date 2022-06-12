% A comparison of the Ridgecrest and Hector Mine ruptures 
clear all;
close all;
fig1 = figure(1);
clf

fault_tr=1; % put fault_tr=1 if want to plot the geologic fault trace
pr=0;       % put pr=1 if want to print output into files
units=0;

bgrey = [1 1 1]*0.7;
xlim=[-25 20]; ylim=[-3 58];
Fontsize=12;

%Ridgecrest

xo=-117.5; yo=35.5;  % origin
[X,Y]=utm2ll(xo,yo,0,1);
xo=X; yo=Y;
x1=-117.599; y1=35.770; % epicenter
[X,Y]=utm2ll(x1,y1,0,1);
x1=(X-xo)*1e-3;
y1=(Y-yo)*1e-3;

if fault_tr==1
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
% [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');
 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
% [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'m');
 [i]=plot_faults(faults,dim,xlim,ylim,'m');
 plot(x1,y1,'rp','MarkerSize',14,'MarkerFaceColor','r','MarkerEdgeColor','k'), hold on  % epicenter location

end

load foc.mech;
M=foc(:,11);
b=find(M>2);
lat=foc(b,8);
lon=foc(b,9);

[X,Y]=utm2ll(lon,lat,0,1);
x1=(X-xo)*1e-3;
y1=(Y-yo)*1e-3;

plot(x1,y1,'bs','MarkerSize',12,'MarkerFaceColor','b','MarkerEdgeColor','c'), hold on  % epicenter location

%   text(-15,2,'M7.1 Ridgecrest','Fontsize',14,'Color','k'), hold on

axis('equal')
xlabel('Eastings, km','Fontsize',14);
ylabel('Northings, km','Fontsize',14);
set(gca,'Fontsize',14,'XLim',xlim,'YLim',ylim,'box','on');

orientation = [0.05, 0.1, 7, 7; 0.1, 0.1,7, 7];
computerNames = {'foxbat'; 'flagon'};
[~, thiscomputer] = system('hostname');
option = strcmp(deblank(thiscomputer), computerNames);
%set(gcf,'defaulttextinterpreter','none','Visible','off','PaperOrientation', 'landscape','PaperPositionMode','manual', 'PaperPosition',orientation(option, :))


if pr==1
 print('-dpng','-r300','RC_focmech')
end

