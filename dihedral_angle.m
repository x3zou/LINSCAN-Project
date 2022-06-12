%Calculate the dihedral angles
clear 
clc

adjust2 = 180; % degrees
bin=5; % bin size for strikes from the slip model , deg
Fs=12; % font size

% read model geometry (segment coords/moment tensor)
load date1dihedral_input.txt;
x=date1dihedral_input;

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
er=[]; % error of difference
ij=[];
k=0;
figure(4)
clf
for i=1:length(s)
 for j=i:length(s)
%  if i~=j & ((s(i)>90 & s(j)<90) | (s(i)<90 & s(j)>90))
  if i~=j & (abs(abs(r3(i))-abs(r3(j)))>90)

      dff=abs(s4(i)-s4(j)); 
      df=[df dff];
      er=[er err(i)+err(j)];
%if dff >95, fprintf('%d %d \n',i,j); end
%      if dff >95
      if dff <31
        k=k+1; 
        ij(k,1)=i;  
        ij(k,2)=j;  
      end

  end
 end
end
%ind=find(df>180);
%df(ind)=df(ind)-180;

%cmp=unique(ij);
cmp=ij;
save flt_num cmp -ascii;

x=0:bin:120; % angle bins

del=0.25; % increment for numerical integration

%for j=1:length(x)-1
% bn=find(df>x(j) & df<=x(j+1));
% sm(j)=length(bn);
% ers(j)=std(df(bn));
%end


sm_err=0*x;
qm=0;
for j=1:length(x)-1
 pij=[];
 cnt1=0;
 cnt2=0;
 bn=find(df>x(j) & df<=x(j+1));
 sm(j)=length(bn);
% ers(i)=std(er(bn))/sqrt(length(bn));
% ers(j)=std(df(bn));
% ers(i)=3*std(df(bn))/sqrt(length(bn));
% bne=find(df>x(i)-ers(i) & df<=x(i+1)+ers(i));
 for i=1:sm(j)
%  [pij(i)]=gauss_int(df(bn(i)),er(bn(i)),x(1):del:x(length(x)));
  [pij(i)]=gauss_int(df(bn(i)),er(bn(i)),x(j):del:x(j+1));
  if sm(j)>0
   qi=pij(i)*del/sqrt(2*pi);
   qm=max([qm qi]);
   cnt1=cnt1+qi;  % Expectation
   cnt2=cnt2+qi*(1-qi);  % Bernoulli variance
  end
 end 

fprintf('Sum of probabilities: %f \n',cnt1);
fprintf('max probability: %f \n',qm);
% sm(j)=cnt1;
 sm_err(j)=sqrt(cnt2)/cnt1*sm(j);  % standard deviation
fprintf('Error: %f , bin: %d \n',sm_err(j),j); 
%   sm_err(j)=sm(j)*cnt2/cnt1;
end

   max(sm_err(:))

%return

clr=[0.7 0.7 0.7];
figure(5)
set(gcf, 'PaperPosition',[0.1 0.1 5 4]);
clf
colororder({'k','b'})
%yyaxis left

tik=0.7;
for i=1:length(x)-1
  verty=[0 sm(i) sm(i) 0];
  vertx=[x(i) x(i) x(i+1) x(i+1)];
  fill(vertx,verty,clr,'EdgeColor','k','lineStyle','-','Marker','none'), hold on
  xi=mean([x(i) x(i+1)]);
  yi1=sm(i)-sm_err(i)/2;
  yi2=sm(i)+sm_err(i)/2;
  line([xi xi],[yi1 yi2],'Color','k','lineStyle','-','LineWidth',1), hold on
  line([xi-tik xi+tik],[yi1 yi1],'Color','k','lineStyle','-','LineWidth',1), hold on
  line([xi-tik xi+tik],[yi2 yi2],'Color','k','lineStyle','-','LineWidth',1), hold on
end

xlim=[20 110];
yl=get(gca,'ylim');
%xi=30;
%line([xi xi],yl,'Color','r','lineStyle','--','LineWidth',1), hold on
%xi=100;
%line([xi xi],yl,'Color','r','lineStyle','--','LineWidth',1), hold on
%text(32, 100,'$2\theta_1$','FontSize', 16, 'Interpreter','latex'), hold on
%text(90, 100,'$2\theta_2$','FontSize', 16, 'Interpreter','latex'), hold on

%histogram(df);
xlabel('Dihedral angle, deg.');
ylabel('Number of conjugate pairs');

%yyaxis right
%x1=20:90;
%line(x1,1./tand(x1),'lineStyle','-','lineWidth',2,'Color','b'), hold on
%ylabel('Coefficient of friction \mu')     
%hold off
%ax=gca;
%ax.YAxis.Color = 'r';
%set(gca,'xlim',xlim,'ylim',[0 1],'box','on','FontSize',Fs), hold on

set(gca,'FontSize',Fs,'box','on');

print(gcf,'-dpng','-r300','hist');

function [p] = gauss_int(xi,si,z)
 % Compute the integral
  p=0;
for i=1:length(z)
  p=p+exp(-(xi-z(i))^2/2/si^2)/si;
end

end
