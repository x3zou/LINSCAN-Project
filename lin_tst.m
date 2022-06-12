clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
% %minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
% %minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
Fs=14;

deg2rad=pi/180.;


pr=1;
fault_tr = 1; % flag for plotting the fault trace
% regional
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);
 grey = [1 1 1]*0.6;

D = 0.5; % distance threshold for selecting nearby points
style = ':';
color = 'b';
width = 2;
color = [rand(1) rand(1) rand(1)];
X=[]; Y=[]; Len=[]; S=[];

x=[-49.3247  27.5632 -50.2335  25.1130];
s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
 ind=find(s>0 & s<90);
 s(ind)=90-s(ind);
 ind=find(s>90);
 s(ind)=450-s(ind);
 ind=find(s<0 & s>-90);
 s(ind)=270-s(ind);
 ind=find(s<-90);
 s(ind)=-90-s(ind);

s

x=[-50.2335  25.1130 -49.3247  27.5632];
s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
 ind=find(s>0 & s<90);
 s(ind)=90-s(ind);
 ind=find(s>90);
 s(ind)=450-s(ind);
 ind=find(s<0 & s>-90);
 s(ind)=270-s(ind);
 ind=find(s<-90);
 s(ind)=-90-s(ind);

s



return

f = 'lin_segments.dat';
%f = 'lin_segments.rev';
if isfile(f)
 fid = fopen(f,'r');
 load(f);
 [m,n]=size(lin_segments);
 x=lin_segments;
 X=0.5*(x(:,1)+x(:,3));
 Y=0.5*(x(:,2)+x(:,4));
 Len=sqrt((x(:,4)-x(:,2)).^2+(x(:,3)-x(:,1)).^2);
 s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
 
 ind=find(s>0 & s<90);
 s(ind)=90-s(ind);
 ind=find(s>90);
 s(ind)=450-s(ind);
 ind=find(s<0);
 s(ind)=270-s(ind);
 ind=find(s>360);
 s(ind)=450-s(ind);
 S=s;

% for i=1:m
%  line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
% end
 fclose(fid);
end

f = 'fc.dat';
if isfile(f)
 fid = fopen(f,'r');
 load(f);
 x=fc;
 [m,n]=size(x);
end
 z=x;


 for i=1:m
   figure(1)
clf
  line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
  line([z(i,1) z(i,3)],[z(i,2) z(i,4)],'lineStyle','-','Color','r','LineWidth',1), hold on
pause
 end

return

%pause

fid=fopen('fc.dat','wt');

% % plot all of the clusters:
x0=X; y0=Y; L=Len; strike=S;
%for i=1:0
for i=1:length(x0)
%for i=length(x0):-1:1
 color = [rand(1) rand(1) rand(1)];  %  fprintf(fid, '%s \n',color);

 color = [0 0 0];  
 %figure(1)
 %clf
 plot_fm(i,x0(i),y0(i),x(i,:),strike(i),L(i),D,color,style,width,fid);
end

fclose(fid);


axis('equal');
xlabel('Eastings, km');
ylabel('Nortings, km');
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);

saveas(gcf, 'lin_fm', 'png')
return


function [x1,y1,x2,y2] = plot_line(x,y,z,strike,L,color,style,width,fid)
 % Compute end points
  x1 = x + L/2*sind(strike);
  x2 = x - L/2*sind(strike);
  y1 = y + L/2*cosd(strike);
  y2 = y - L/2*cosd(strike);
  if fid > 0
   fprintf(fid,'%f %f \n',x1,y1);
   fprintf(fid,'%f %f \n',x2,y2);
   fprintf(fid,'> \n',color);
  end
%  fprintf(fid, '%s \n',color);
  fprintf('%f \n',strike);

 for i=1:length(z(:,1))
  line([z(i,1) z(i,3)],[z(i,2) z(i,4)],'lineStyle','-','Color','b','LineWidth',2), hold on
 end

  line([x1 x2],[y1 y2],'lineWidth',1,'lineStyle','-','Color','r'), hold on

end




function plot_fm(k,x0,y0,xlin,strike,L,D,color,style,width,fid)
    figure(1)
clf
%[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,fid);
[x1,y1,x2,y2]=plot_line(x0,y0,xlin,strike,L,color,style,width,0);

plot([x1 x2],[y1 y2],'ko'), hold on
%pause
end

