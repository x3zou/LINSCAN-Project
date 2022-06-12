clear all;
close all;

Fs=14;
fault_tr=1;


x1=[];
load /Users/fialko/w/proj/SCal_eq_2019/qtm_final_12dev.dat;
x=qtm_final_12dev;
figure
X=x(:,9); % lon
Y=x(:,8); % lat
Z=x(:,10);% depth
%plot(X, Y,'.k'), hold on
%plot3(X, Y, -Z,'.k'), hold on
%return

%yr=x(:,1); % year
%d=x(:,3);  % day
%h=x(:,4);  % hour
%m=x(:,5);  % min

event = find(X> -118.2 & X < -117 & Y < 36.3 & Y > 35.3);

orig=[-117.5 35.5];
x0=orig(1); y0=orig(2);
[xo,yo]=utm2ll(x0,y0,0,1);

%plot(X(event), Y(event),'.k'), hold on
%plot(X(event), Y(event),'.k'), hold on
%plot3(X(event), Y(event), -Z(event),'.k'), hold on


%plot(s(:,1),s(:,2),'or'), hold on
%plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on


%saveas(gcf,'seism.fig')
%return

x1(:,1) =X(event);  
x1(:,2) =Y(event);
x1(:,3) =-Z(event);


%plot3(x1(:,1), x1(:,2), x1(:,3),'.k'), hold on


if fault_tr==1
 [x1(:,1),x1(:,2)]=utm2ll(x1(:,1),x1(:,2),0,1);
 x1(:,1)=(x1(:,1)-xo)*1e-3;
 x1(:,2)=(x1(:,2)-yo)*1e-3;
end

xmin=min(x1(:,1)); xmax=max(x1(:,1)); ymin=min(x1(:,2)); ymax=max(x1(:,2));
xl=[xmin xmax]; yl=[ymin ymax];


plot3(x1(:,1), x1(:,2), x1(:,3),'.k'), hold on

save preeq_ross.utm x1 -ascii
return

col='m';
if fault_tr==1
  Hector_fault_ll=[];  Landers_fault_ll=[]; ca_faults_ll=[];
  Hector_dim=[]; Landers_dim=[]; ca_dim_ll=[];
%  load /Users/fialko/w/dat/Volumes/scratch/work/insar/landers/pic/Landers_fault_ll.dat;
%  load /Volumes/scratch/work/insar/landers/pic/Landers_dim.dat;
%  load /Volumes/scratch/work/insar/landers/pic/Hector_dim.dat;
%  load /Volumes/scratch/work/insar/landers/pic/Hector_fault_ll.dat;
  load /Users/fialko/w/dat/ca_dim_ll.dat;
  load /Users/fialko/w/dat/ca_faults_ll.dat;
  faults=[Hector_fault_ll' Landers_fault_ll' ca_faults_ll']';
  dim=[Hector_dim' Landers_dim' ca_dim_ll']';
%  faults=[Hector_fault_ll']';
%  dim=[Hector_dim']';
 [sf,l]=size(dim);
% [sfhm,lhm]=size(Hector_dim);
 [sfhm,lhm]=size(dim);
end

if fault_tr==1
 [faults(:,1),faults(:,2)]=utm2ll_kludge(faults(:,1),faults(:,2),x0,y0);
 [i]=plot_faults(faults,dim,xl,yl,col);
end %faults

x=[];
% plot a sphere at a given point
x(1,1)=-116.087; x(1,2)=33.663; 
[x1,y1]=utm2ll_kludge(x(:,1),x(:,2),x0,y0);
[x,y,z] = sphere(100); hold on
R=1; % sphere radius, km
hs1 = surf(x*R+x1,y*R+y1,z*R+0,'FaceAlpha',0.9); shading flat
q1 = get(hs1);
set(hs1, 'FaceColor', [0 1 1]); 
hold off

x=[];
% plot a sphere at a given point
x(1,1)=-116.0213; x(1,2)=33.6068; 
[x1,y1]=utm2ll_kludge(x(:,1),x(:,2),x0,y0);
[x,y,z] = sphere(100); hold on
%R=0.5; % sphere radius, km
hs2 = surf(x*R+x1,y*R+y1,z*R+0,'FaceAlpha',0.9); shading flat
q2 = get(hs2);
set(hs2, 'FaceColor', [1 0 0]);

x=[];
% plot a sphere at a given point
x(1,1)=-115.93; x(1,2)=33.533; 
[x1,y1]=utm2ll_kludge(x(:,1),x(:,2),x0,y0);
[x,y,z] = sphere(100); hold on
%R=0.5; % sphere radius, km
hs3 = surf(x*R+x1,y*R+y1,z*R+0,'FaceAlpha',0.9); shading flat
q3 = get(hs3);
set(hs3, 'FaceColor', [0 0 1]);


set(gca,'box','on'), hold on
grid on, hold on
set(gca,'FontSize',Fs), hold on
axis equal;

saveas(gcf,'micro_eq.png');

