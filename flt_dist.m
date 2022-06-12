clear all
close all

adjust = 90; % degrees
adjust2 = 180; % degrees
bin=5; % bin size for strikes from the slip model , deg
Fs=14; % font size

% read model geometry (segment coords/moment tensor)
load fc2.new;
x=fc2;
ind=find((x(:,4)+x(:,2))/2 > -40);
x=x(ind,:);

s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
Len=sqrt((x(:,4)-x(:,2)).^2+(x(:,3)-x(:,1)).^2);
s4=s;
 ind=find(s>0 & s<90);
 s(ind)=90-s(ind);
 ind=find(s>90);
 s(ind)=450-s(ind);
 ind=find(s<0 & s>-90);
 s(ind)=270-s(ind);
 ind=find(s<-90);
 s(ind)=-90-s(ind);
st=s;

% 1st nodal plane:
s1=x(:,5);
d1=x(:,6);
r1=x(:,7);
% 2nd nodal plane:
s2=x(:,8);
d2=x(:,9);
r2=x(:,10);

err=x(:,11);

ds=abs(st-s1);
s3=0*s1; d3=0*d1; r3=0*r1; 
for i=1:length(ds)
if abs(rem(round(ds(i)/90),2))==0
  s3(i)=s1(i);
  d3(i)=d1(i);
  r3(i)=r1(i);
else
  s3(i)=s2(i);
  d3(i)=d2(i);
  r3(i)=r2(i);
end
end

deg2rad=pi/180.;
 figure(2)
% plot(st,x(:,5),'ko'), hold on
 plot(1:length(ds),ds,'ko'), hold on

% for i=1:length(x(:,1)) 
% figure(4)
% clf
% line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','lineWidth',2,'Color','r'), hold on
% text(mean([x(i,1) x(i,3)])-2,mean([x(i,2) x(i,4)])-2, ['strike: ' num2str(s3(i))]), hold on
% text(mean([x(i,1) x(i,3)])-2,mean([x(i,2) x(i,4)])-2.5, ['rake: ' num2str(r3(i))]), hold on
% text(mean([x(i,1) x(i,3)])-2,mean([x(i,2) x(i,4)])-3, ['dif: ' num2str(ds(i))]), hold on
% %  bb([s3(i) d3(i) r3(i)], mean([x(i,1) x(i,3)])-4,mean([x(i,2) x(i,4)])-4, 1.5, 0, 'r'), hold on
%   [fm]=sdr2mij([s3(i) d3(i) r3(i)],deg2rad);
% focalmech(fm, mean([x(i,1) x(i,3)])-4,mean([x(i,2) x(i,4)])-4, 3, 'r','dc'), hold on
% axis equal;
%   pause         
% end
% 
% return

figure(3)
%fig


for i=1:length(x(:,1))
line([x(i,1) x(i,3)],[x(i,2) x(i,4)],'lineStyle','-','lineWidth',2,'Color','r'), hold on
end
axis equal;
%return

k=[];
df=[]; % difference
ds=[]; % difference
er=[]; % error of difference
figure(4)
clf
for i=1:length(s)
 for j=i:length(s)
%  if i~=j & ((s(i)>90 & s(j)<90) | (s(i)<90 & s(j)>90))
  if i~=j & (abs(abs(r3(i))-abs(r3(j)))<90)

      dff=abs(s4(i)-s4(j));
      dst=sqrt((x(i,4)+x(i,2)-x(j,4)-x(j,2))^2+(x(i,3)+x(i,1)-x(j,3)-x(j,1))^2)/2;
      df=[df dff];
      ds=[ds dst];
  end
 end
end
%ind=find(df>180);
%df(ind)=df(ind)-180;


x=0:5:50;
cnt=0*x;
sm=0*x;

for i=1:length(x)-1
 bn=find(df>x(i) & df<=x(i+1));
 if length(bn)>0
  sm(i)=min(ds(bn));
  cnt(i)=length(bn);
 end
end

clr=[0.7 0.7 0.7];
figure(5)
set(gcf, 'PaperPosition',[0.1 0.1 5 4]);
clf
%colororder({'k','r'})
%yyaxis left

plot(ds,df,'k.','MarkerSize',10), hold on
xlabel('Distance between faults, km');
ylabel('Difference in fault strikes, deg');
set(gca,'xlim',[0 80],'ylim',[0 45],'FontSize',Fs), hold on
print(gcf,'-dpng','-r300','flt_dist');

return

