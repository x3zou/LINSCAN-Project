clear all;
close all;
fault_tr = 1; % flag for plotting the fault trace
% regional
minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
Fs=14;
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
grey = [1 1 1]*0.5;
%grey = [0.4,0.4,0.4];

load ps_off.dat;
x=ps_off;

figure(1)

if fault_tr==1
  Ridgecrest_fault_ll=[];  ca_faults_ll=[];
  Ridgecrest_dim=[];  ca_dim_ll=[];
  load Ridgecrest7_fault_ll.dat;
  load Ridgecrest7_dim.dat;
  faults=Ridgecrest7_fault_ll;
  dim=Ridgecrest7_dim;
  [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
  faults(:,1)=(faults(:,1)-xo)*1e-3;
  faults(:,2)=(faults(:,2)-yo)*1e-3;
% [sf,l]=size(dim);
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'r');

 load Ridgecrest2.trace;
 faults=Ridgecrest2;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3; 
 line(faults(:,1),faults(:,2),'lineWidth',3,'lineStyle','-','Color','b'), hold on

end

[x(:,1),x(:,2)]=utm2ll(x(:,1),x(:,2),0,1);
x(:,1)=(x(:,1)-xo)*1e-3;
x(:,2)=(x(:,2)-yo)*1e-3;

quiver(x(:,1),x(:,2),x(:,4),x(:,3),0), hold on

f = 'gps.dat';
T=readtable([f],'Delimiter',' ','ReadVariableNames',0);
C=table2cell(T);
wid=length(C(1,:));
N=cell2mat(C(:,1:2)); % coordinates
L=C(:,3);
x=N(:,1); y=N(:,2);
[xx,yy]=utm2ll(x(:),y(:),0,1);
xu=reshape((xx-xo)*1e-3,size(x));
yu=reshape((yy-yo)*1e-3,size(x));
for i=1:length(xu(:))
 xsite=char(L(i));
 site=xsite(2:5);
 text(xu(i),yu(i),site), hold on
end

% text(-45,5,'1984-2011','FontSize',Fs);

axis('equal');
xlabel('Eastings, km');
ylabel('Nortings, km');
%set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
%legend('earthquakes','\sigma_1','$\dot{\epsilon}$', 'Interpreter','latex');

saveas(gcf, 'ps_disp', 'png')
return
