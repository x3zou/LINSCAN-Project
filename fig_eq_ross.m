clear all;
close all;

load ridgecrest_qtm.cat;
x=ridgecrest_qtm;
figure
X=x(:,9);
Y=x(:,8);
Z=x(:,10);
plot(X, Y,'.k','MarkerSize',1), hold on

% figure out the correct lat/lon aspect ratio
del=0.01;
orx=mean(X(:)); ory=mean(Y(:));
[tx1,ty1]=utm2ll(orx,ory,0,1);
[tx2,ty2]=utm2ll(orx+del,ory,0,1);
[tx3,ty3]=utm2ll(orx,ory+del,0,1);
llrat=[1 abs((tx2-tx1)/(ty3-ty1)) 1];

%load eqn2.ross;
%x=eqn2;
%plot(x(:,1), x(:,2),'.c'), hold on

load Ridgecrest1.trace;
s=Ridgecrest1;
%load Ridgecrest2.trace;
%s=Ridgecrest2;

plot(s(:,1),s(:,2),'.g','MarkerSize',10), hold on

load Ridgecrest2.trace;
s=Ridgecrest2;
plot(s(:,1),s(:,2),'.g','MarkerSize',10), hold on

%plot(s(:,1),s(:,2),'ob'), hold on
%plot3(s(:,1),s(:,2),0*s(:,1),'or'), hold on
ylabel('Longtitude, deg'), hold on
xlabel('Latitude, deg.'), hold on
set(gca,'xlim',[-117.9 -117.3],'ylim',[35.5 36.1]), hold on
set(gca,'box','on','DataAspectRatio',llrat), hold on
set(gca,'FontSize',14), hold on

%saveas(gcf,'seism_ross.fig')
saveas(gcf,'seism_ross.png')
%return
