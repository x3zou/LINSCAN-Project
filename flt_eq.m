clear all;
close all;

load eq1995.dat;
x=eq1995;
X=x(:,1);
Y=x(:,2);
Z=x(:,3);
%plot(X, Y,'.k'), hold on
orig=[-117.5 35.5];
xo=orig(1); yo=orig(2);
[xo,yo]=utm2ll(xo,yo,0,1)

event = find(X> -117.8 & X < -117.5 & Y < 35.85 & Y > 35.65);
x=X(event); 
y=Y(event);
[X1,Y1]=utm2ll(x,y,0,1);
xU=(X1-xo)*1e-3;
yU=(Y1-yo)*1e-3;
 

plot3(xU, yU, -Z(event),'.k'), hold on
axis('equal'), hold on

return

load Hauksson_ridgecrest_17_July.reloc;
x=Hauksson_ridgecrest_17_July;
figure
X=x(:,3);
Y=x(:,2);
Z=x(:,4);
%plot(X, Y,'.k'), hold on

d=x(:,13);
h=x(:,14);
m=x(:,15);

  event = find(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84));
%event = find(~(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84)));

%event = find(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84) & d < 5 & (d == 5 & h < 20) );

%plot(X(event), Y(event),'.k'), hold on
%plot(X(event), Y(event),'.k'), hold on
plot3(X(event), Y(event), -Z(event),'.k'), hold on

load Ridgecrest1.trace;
s=Ridgecrest1;

%plot(s(:,1),s(:,2),'or'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

load Ridgecrest2.trace;
s=Ridgecrest2;

%plot(s(:,1),s(:,2),'ob'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

saveas(gcf,'seism.fig')
return

x1(:,1) =X(event);  
x1(:,2) =Y(event);

%save eq2.dat x1 -ascii
