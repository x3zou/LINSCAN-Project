clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
% %minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
% %minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
grey = [1 1 1]*0.6;
grey = [1 1 1]*0;
Fs=14;

load plane.out;
load fmec.dat;

figure(1)
clf

load eq_gc.utm;
x=eq_gc;
plot3(x(:,1),x(:,2),x(:,3),'LineStyle','none','Marker','.','MarkerSize',1,...
    'MarkerEdgeColor','m','MarkerFaceColor','m'), hold on;

%plot3(fmec(:,1),fmec(:,2),fmec(:,3),'LineStyle','none','Marker','.','MarkerSize',1,...
%    'MarkerEdgeColor',grey,'MarkerFaceColor',grey), hold on
%plot3(plane(:,4),plane(:,5),plane(:,6),'LineStyle','none','Marker','o','MarkerSize',5,...
%    'MarkerEdgeColor','b','MarkerFaceColor','b'), hold on


f = 'lin_segments.dat2';
fid = fopen(f,'wt');

%for i=1:0
for i=1:length(plane(:,1))
  xb= [plane(i,4) plane(i,1) plane(i,7) plane(i,10)];    
  yb= [plane(i,5) plane(i,2) plane(i,8) plane(i,11)];    
  zb= [plane(i,6) plane(i,3) plane(i,9) plane(i,12)];    
%  fill3(xb,yb,zb,1,'EdgeColor','k','LineWidth',1,'FaceAlpha',0.5,'EdgeAlpha',1); hold on  

  x=[mean([plane(i,4) plane(i,1)]) mean([plane(i,7) plane(i,10)])]; 
  y=[mean([plane(i,5) plane(i,2)]) mean([plane(i,8) plane(i,11)])];
  line([x(1) x(2)],[y(1) y(2)],[0 0],'LineWidth',3,'Color','k'), hold on 

  fprintf(fid,'%9.4f %9.4f %9.4f %9.4f\n',x(1), y(1), x(2), y(2));
  
end

load plane.out2;
for i=1:length(plane(:,1))
  xb= [plane(i,4) plane(i,1) plane(i,7) plane(i,10)];    
  yb= [plane(i,5) plane(i,2) plane(i,8) plane(i,11)];    
  zb= [plane(i,6) plane(i,3) plane(i,9) plane(i,12)];    
%  fill3(xb,yb,zb,2,'EdgeColor','k','LineWidth',1,'FaceAlpha',0.5,'EdgeAlpha',1); hold on
  x=[mean([plane(i,4) plane(i,1)]) mean([plane(i,7) plane(i,10)])]; 
  y=[mean([plane(i,5) plane(i,2)]) mean([plane(i,8) plane(i,11)])];
  line([x(1) x(2)],[y(1) y(2)],[0 0],'LineWidth',3,'Color','k'), hold on 

  fprintf(fid,'%9.4f %9.4f %9.4f %9.4f\n',x(1), y(1), x(2), y(2));
end
fclose(fid);


axis equal;
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'zlim',[-12 0],'box','on','FontSize',Fs);
%text(-1,3,num2str(s)), hold on

return


