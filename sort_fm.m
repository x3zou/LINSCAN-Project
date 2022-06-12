clear all;
close all;

load foc.mech;
%M=foc(:,11);
%b=find(M>6);
Y=foc(:,8);
X=foc(:,9);
Z=foc(:,10);
d=foc(:,3);
h=foc(:,4);
m=foc(:,2);
mn=foc(:,5);
M=foc(:,11);

t=m/12+d/365+h/365/24+mn/365/24/60;
tm6=7/12+4/365+17/365/24+33/365/24/60;
tm7=7/12+6/365+3/365/24+19/365/24/60;

%plot(X, Y,'.k'), hold on

orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1);

%event = find(X> -117.8 & X < -117.1 & Y < 36 & Y > 35.5);
event = find(X> -117.8 & X < -117.1 & Y < 36 & Y > 35.5 & ~(X>-117.56 & Y > 35.84) & t > tm6 & t < tm7 & M > 0);
M=M(event);
T=t(event);
x=X(event); 
y=Y(event);
[X1,Y1]=utm2ll(x,y,0,1);
xU=(X1-xo)*1e-3;
yU=(Y1-yo)*1e-3;
 
%plot3(xU, yU, -Z(event),'.k'), hold on
for i=1:length(event)
  plot3(xU(i), yU(i), -Z(event(i)),'.k','MarkerSize',M(i)*6), hold on
end

  plot3(-3, 25, -7,'.r','MarkerSize',50), hold on

axis('equal'), hold on
xlabel('Eastings, km','Fontsize',14);
ylabel('Northings, km','Fontsize',14);
box on
grid on

return

%Omori's law

t0=min(T); t1=max(T);
nt=20;
bins=t0:(t1-t0)/nt:t1;

for i=1:length(bins)-1
 b=find(T>bins(i) & T < bins(i+1));
 numq=length(b); 
end

figure(2)
plot(bins,numq,'.k'), hold on

return


load Ridgecrest1.trace;
s=Ridgecrest1;

%plot(s(:,1),s(:,2),'or'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

load Ridgecrest2.trace;
s=Ridgecrest2;

%plot(s(:,1),s(:,2),'ob'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

%saveas(gcf,'seism.fig')
return

x1(:,1) =X(event);  
x1(:,2) =Y(event);

%save eq2.dat x1 -ascii
