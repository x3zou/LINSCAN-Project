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

load lin.utm;
xlin=lin;

figure(1)
line(xlin(:,1),xlin(:,2),'lineStyle','none','Marker','.','MarkerSize',0.1,'MarkerFaceColor',grey,'MarkerEdgeColor',grey), hold on
%return

f = 'gpss.dat';
T=readtable([f],'Delimiter',' ','ReadVariableNames',0);
C=table2cell(T);
wid=length(C(1,:));
N=cell2mat(C(:,1:2)); % coordinates
L=C(:,3);             % station labels
x=N(:,1);
y=N(:,2);
[xx,yy]=utm2ll(x(:),y(:),0,1);
xu=reshape((xx-xo)*1e-3,size(x));
yu=reshape((yy-yo)*1e-3,size(x));
line(xu,yu,'lineStyle','none','Marker','^','MarkerSize',10,'MarkerFaceColor','c','MarkerEdgeColor','k'), hold on
text(xu+0.3,yu-1.5,char(L)), hold on

%fname = 'gnss.utm';
%%arr = [xu yu char(L)];
%arr = [xu yu];
%save(fname,'arr','-ascii');
%return

f = 'coso.dat';
T=readtable([f],'Delimiter',' ','ReadVariableNames',0);
C=table2cell(T);
wid=length(C(1,:));
N=cell2mat(C(:,1:2)); % coordinates
L=C(:,3);             % station labels
x=N(:,1);
y=N(:,2);
[xx,yy]=utm2ll(x(:),y(:),0,1);
xu=reshape((xx-xo)*1e-3,size(x));
yu=reshape((yy-yo)*1e-3,size(x));
line(xu,yu,'lineStyle','none','Marker','^','MarkerSize',10,'MarkerFaceColor','y','MarkerEdgeColor','k'), hold on
text(xu-3.5,yu-1.5,char(L)), hold on
%fname = 'coso.utm';
%%arr = [xu yu char(L)];
%arr = [xu yu];
%save(fname,'arr','-ascii');

%return

% principal stress axis:
style = '-';
color = 'k';
width = 1;
L = 5;
load hmax_G10N30_rect_gmt.dat;
x=hmax_G10N30_rect_gmt;
[x(:,1),x(:,2)]=utm2ll(x(:,1),x(:,2),0,1);
x(:,1)=(x(:,1)-xo)*1e-3;
x(:,2)=(x(:,2)-yo)*1e-3;
%for i=1:length(x(:,1))
%  plot_line(x(i,1),x(i,2),x(i,5),L,color,style,width);
%end

del=2;
[X,Y]=meshgrid(min(x(:,1)):del:max(x(:,1)),min(x(:,2)):del:max(x(:,2)));
F = scatteredInterpolant(x(:,1),x(:,2),x(:,5),'linear');
Z = F(X,Y);
%interp2(x(:,1),x(:,2),x(:,5),X(:),Y(:)); 
del=5;
sx = minlon-20:del:maxlon;
sy=-5*ones(size(sx));
%sy=(1:length(sx))*(-20); %*minlat;
h=streamline(stream2(X,Y,sind(Z),cosd(Z),sx,sy));
% save in a file
fid1 = fopen('stress_lines.utm','a');
fmt = ' %10.7f %10.7f \n';
for i = 1:length(h)
 for j = 1:length(h(i).XData)
  fprintf(fid1,fmt, [h(i).XData(j), h(i).YData(j)]);
 end
 fprintf(fid1,'> \n');
end
fclose(fid1);

set(h,'Color','k','LineWidth',0.5,'LineStyle','--');
ind = find(X > -15 & X < -5 & Y > 20 & Y < 40);
mean(Z(ind))
% principal axis of strain rate:
load p_strain.dat;
x=p_strain;
F = scatteredInterpolant(x(:,1),x(:,2),x(:,3),'linear');
Z = F(X,Y);
h=streamline(stream2(X,Y,sind(Z),cosd(Z),sx,sy));

fid1 = fopen('strain_lines.utm','a');
fmt = ' %10.7f %10.7f \n';
for i = 1:length(h)
 for j = 1:length(h(i).XData)
  fprintf(fid1,fmt, [h(i).XData(j), h(i).YData(j)]);
 end
 fprintf(fid1,'> \n');
end
fclose(fid1);

set(h,'Color','m','LineWidth',0.5,'LineStyle','-');
%ind = find(X > -20 & X < 10 & Y > 10 & Y < 40);
ind = find(X > -15 & X < -5 & Y > 20 & Y < 40);
mean(Z(ind))

strike = 42;
style = ':';
color = 'b';
width = 2;

fid=fopen('LL.dat','wt');

D = 2; % distance threshold for selecting nearby points

% 1st cluster
x0 =-50; y0 = 25;
L = 20;
[x1,y1,x2,y2]=plot_line(x0,y0,strike,L,color,style,width,fid);
plot([x1 x2],[y1 y2],'ko'), hold on
v1=[x1 y1 0];  v2=[x2 y2 0]; 

nx=[]; ny=[];
pts=[xlin(:,1) xlin(:,2) 0*xlin(:,1)];

[ind] = point_to_line(pts, v1, v2, D);
%ind=find(nearby<D);
%ind=find(~isnan(nearby));
nx = xlin(ind,1);
ny = xlin(ind,2);

plot(nx,ny,'r.'), hold on

return





x0 =-11; y0 = 19;
L = 9;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-4; y0 = 35;
L = 7;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-7; y0 = 60;
L = 9;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-3.5; y0 = 58;
L = 9;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-12.751225; y0 = 62.837010;
L = 17;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-36.158088; y0 = 69.699755;
L = 5;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-31; y0 = 63;
L = 15;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-31; y0 = 63;
L = 15;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =17; y0 = 52.5;
L = 6;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-5.5; y0 = 15;
L = 7;
plot_line(x0,y0,strike,L,color,style,width,fid);

fclose(fid);

fid=fopen('RL.dat','wt');

strike = -40;
style = ':';
color = 'r';

x0=-40; y0= 48.28;
L = 10;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0=-12.230392; y0= 27.787990;
L = 6;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0=10.379902; y0=53.768382;
L = 7;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-14; y0 = 51;
L = 10;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-25.5; y0 = 55;
L = 10;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-2; y0 = 60;
L = 5;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =19; y0 = 50;
L = 5;
plot_line(x0,y0,strike,L,color,style,width,fid);

x0 =-47; y0 = 31;
L = 5;
plot_line(x0,y0,strike,L,color,style,width,fid);

fclose(fid);

if fault_tr==1
%  Ridgecrest_fault_ll=[];  ca_faults_ll=[];
%  Ridgecrest_dim=[];  ca_dim_ll=[];
%  load Ridgecrest7_fault_ll.dat;
%  load Ridgecrest7_dim.dat;
%  faults=Ridgecrest7_fault_ll;
%  dim=Ridgecrest7_dim;
%  [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
%  faults(:,1)=(faults(:,1)-xo)*1e-3;
%  faults(:,2)=(faults(:,2)-yo)*1e-3;
%% [sf,l]=size(dim);
% [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'r');
%
% load Ridgecrest2.trace;
% faults=Ridgecrest2;
% [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
% faults(:,1)=(faults(:,1)-xo)*1e-3;
% faults(:,2)=(faults(:,2)-yo)*1e-3; 
% line(faults(:,1),faults(:,2),'lineWidth',3,'lineStyle','-','Color','b'), hold on

 load trace_ver.ll;
 faults=trace_ver;
 load trace_ver.dat;
 dim=trace_ver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');
 load trace_notver.ll;
 faults=trace_notver;
 load trace_notver.dat;
 dim=trace_notver;
 [faults(:,1),faults(:,2)]=utm2ll(faults(:,1),faults(:,2),0,1);
 faults(:,1)=(faults(:,1)-xo)*1e-3;
 faults(:,2)=(faults(:,2)-yo)*1e-3;
 [i]=plot_faults(faults,dim,[minlon maxlon],[minlat maxlat],'g');

end

% text(-45,5,'1984-2011','FontSize',Fs);
load preeq_ross.utm;
x=preeq_ross;
line(x(:,1),x(:,2),'lineStyle','none','Marker','.','MarkerSize',0.1,'MarkerFaceColor','m','MarkerEdgeColor','m'), hold on


axis('equal');
xlabel('Eastings, km');
ylabel('Nortings, km');
set(gca,'xlim',[minlon maxlon],'ylim',[minlat maxlat],'box','on','FontSize',Fs);
%legend('earthquakes','\sigma_1','$\dot{\epsilon}$', 'Interpreter','latex');

saveas(gcf, 'lin_trends', 'png')
return


function [x1,y1,x2,y2] = plot_line(x,y,strike,L,color,style,width,fid)
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


  line([x1 x2],[y1 y2],'lineWidth',width,'lineStyle',style,'Color',color), hold on
end



function [ind] = point_to_line(pt, v1, v2, D)
% function [ind] = point_to_line(pt, v1, v2, D)
% return points within distance D to a finite line segment

 L = norm(v1 - v2);

 v1 = repmat(v1,size(pt,1),1); % 1st end of line point
 v2 = repmat(v2,size(pt,1),1); % 2nd end of line point
 a = v1 - v2;
 b = pt - v2;
 c = pt - v1;
 d1= sqrt(sum(c.*c,2)); % square of distance to one end of the line segment
 d2= sqrt(sum(b.*b,2)); % square of distance to the other end of the line segment

% d = sqrt(sum(cross(a,b,2).^2,2) ./ sum(a.^2,2)); %distance to line
 d = sqrt(sum(cross(a,b,2).^2,2))/L; %distance to line

 ind=find(d < D & d1+d2 < L+2*D); % index to nearby points

end
