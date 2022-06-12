clear all;
close all;

%load eq1995.dat;
%x=eq1995;
%X=x(:,1);
%Y=x(:,2);
%Z=x(:,3);
%%plot(X, Y,'.k'), hold on
%orig=[-117.5 35.5];
%xo=orig(1); yo=orig(2);
%[xo,yo]=utm2ll(xo,yo,0,1)

%event = find(X> -117.8 & X < -117.5 & Y < 35.85 & Y > 35.65);
%x=X(event); 
%y=Y(event);
%[X1,Y1]=utm2ll(x,y,0,1);
%xU=(X1-xo)*1e-3;
%yU=(Y1-yo)*1e-3;
 

%plot3(xU, yU, -Z(event),'.k'), hold on
%axis('equal'), hold on

%return

load ridgecrest_qtm.cat;
x=ridgecrest_qtm;
figure
X=x(:,9);
Y=x(:,8);
Z=x(:,10);
%plot(X, Y,'.k'), hold on
%plot3(X, Y, -Z,'.k'), hold on
%return

d=x(:,3);
h=x(:,4);
m=x(:,5);
mg=x(:,11);
min(mg)

  m=-1:0.5:8;
for i=1:length(m)
  ind=find(mg>m(i));
  num(i)=log10(length(ind));
end
  figure(77)
  plot(m,num,'ko'), hold on

return


%event = find(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84));
%event = find(~(X> -117.8 & Y < 36 & ~(X>-117.56 & Y > 35.84)));

% M5.4
%event = find(X< -117.3 & X> -117.78 & Y < 36 & Y > 35.45 & ~(X>-117.56 & Y > 35.84) & d > 5 | (d == 5 & h > 20) );
event=2084;
X1=X(event); Y1=Y(event); 
event = find((sqrt((X-X1).^2+(Y-Y1).^2) < 0.03)  & ((d == 5 & h > 10) | (d == 6 & h < 3)) );
plot(X(event), Y(event),'sc'), hold on

% M7
%event = find(X< -117.3 & X> -117.78 & Y < 36 & Y > 35.45 & ~(X>-117.56 & Y > 35.84) & d > 5 | (d == 5 & h > 20) );

% M6
%event = find(X< -117.3 & X> -117.78 & Y < 36 & Y > 35.45 & ~(X>-117.56 & Y > 35.84) & d < 5 | (d == 5 & h < 20) );

plot(X(event), Y(event),'.k'), hold on
%event=2084;
%plot(X(event), Y(event),'sc'), hold on
%plot3(X(event), Y(event), -Z(event),'.k'), hold on

load Ridgecrest1.trace;
s=Ridgecrest1;

%plot(s(:,1),s(:,2),'or'), hold on
plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

load Ridgecrest2.trace;
s=Ridgecrest2;

%plot(s(:,1),s(:,2),'ob'), hold on
%plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on

%saveas(gcf,'seism.fig')
%return

x1(:,1) =X(event);  
x1(:,2) =Y(event);
x1(:,3) =-Z(event);

save eq.ross2019M5_4 x1 -ascii
%save eq1.ross x1 -ascii
