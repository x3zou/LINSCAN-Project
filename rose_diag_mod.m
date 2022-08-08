clear
clc

adjust = 90; % degrees
adjust2 = 180; % degrees
R=100; % diagram radius
F=R/3;
bin=5; % bin size for strikes from the slip model , deg
Fs=14; % font size

% read model geometry (strike/moment of the slip patches)

%load lin_segments.dat;
%x=lin_segments;
load ttdihedral_newinput.txt;
x=ttdihedral_newinput;
% ind=find((x(:,4)+x(:,2))/2 >  -40);
% x=x(ind,:);

y1=x(:,2);
y2=x(:,4);
x1=x(:,1);
x2=x(:,3);
for i=1:length(x)
    if y1(i)>y2(i)
        y1(i)=x(i,4);
        y2(i)=x(i,2);
        x1(i)=x(i,3);
        x2(i)=x(i,1);
    end
end
x(:,2)=y1;
x(:,4)=y2;
x(:,1)=x1;
x(:,3)=x2;
s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;

%s=atan2(x(:,4)-x(:,2),x(:,3)-x(:,1))*180/pi;
Len=sqrt((x(:,4)-x(:,2)).^2+(x(:,3)-x(:,1)).^2);
%s=x(:,5);
%r=x(:,7);
s1=Len./Len;
s2=s;

figure(1)
%fig

  ind=find(s<0);
  s(ind)=180+s(ind);

%s=[s Len];
s=[s s1];

% bin the data
%a1=[270:5:360];
a1=[90:5:180];
a2=[0:5:90];

b1=0*a1; b2=0*a2;

for i=1:length(a1)-1
  for j=1:length(s(:,1))
    if s(j,1) >= a1(i) & s(j,1) < a1(i+1)
      b1(i)=b1(i)+s(j,2); 
   end 
  end
end

for i=1:length(a2)-1
  for j=1:length(s(:,1))
    if s(j,1) >= a2(i) & s(j,1) < a2(i+1)
      b2(i)=b2(i)+s(j,2); 
   end 
  end
end

b1m=max(b1); b1=b1/b1m;
b2m=max(b2); b2=b2/b2m;

fprintf('total # of faults: %d \n',sum(b1)*b1m+sum(b2)*b2m);
fprintf('total RL faults: %d ; total LL faults: %d \n',sum(b1)*b1m,sum(b2)*b2m);
fprintf('max RL faults: %d ; max LL faults: %d \n',b1m,b2m);

N=floor(R*1.0); % number of samples for the historgam 

c1=[]; c2=[];
for i=1:length(b1)-1
  av=mean([a1(i) a1(i+1)]);
  for j=1:floor(N*b1(i))
   c1=[c1 av];
  end
end
%c1=-c1+90;

N=floor(R*1); % number of samples for the historgam 

for i=1:length(b2)-1
  av=mean([a2(i) a2(i+1)]);
  for j=1:floor(N*b2(i))
   c2=[c2 av];
  end
end

% fault strikes
sDist = deg2rad(c1);
h2=polarhistogram(sDist, 360/bin, 'BinLimits', [0 2*pi],'FaceColor','r','FaceAlpha',0.3); hold on


sDist = deg2rad(c2);
polarhistogram(sDist, 360/bin, 'BinLimits', [0 2*pi],'FaceColor','b','FaceAlpha',0.3); hold on
ax=gca;
ax.RLim = [0 R];

% 1st set

% M7, strike from focal mech
colorAvgAng = 'r';
sa = 160; % mean angle
sa=-adjust-sa;
%if sa < 0, sa=sa+adjust; end
su = 2;   % upper bound
sl = 2;   % lower bound
lineStyle='-';
% draw_line(sa,colorAvgAng,lineStyle,2);
offset=-0.01;
%draw_arc(sa,sl,su,colorAvgAng,offset,lineStyle);
ax.RLim = [0 R];

% M7, aftershocks
colorAvgAng = 'r';
%sa = -41.55; % mean angle , aftershocks of M6
sa = -39; % mean angle , aftershocks of M6
sa=adjust-sa;
%if sa < 0, sa=sa+adjust; end
su = 0;   % upper bound
sl = 0;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
su = 3;   % upper bound
sl = 3;   % lower bound
lineStyle=':';
% draw_line(sa,colorAvgAng,lineStyle,2);
offset=0.02;
%draw_arc(sa,sl,su,colorAvgAng,offset,lineStyle);
ax.RLim = [0 R];

% M7, strike from moment tensor
colorAvgAng = 'r';
sa = 322; % mean angle
sa=adjust-sa;
%if sa < 0, sa=sa+adjust; end
su = 1;   % upper bound
sl = 1;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
% draw_line(sa,colorAvgAng,'--',2);
offset=-0.01;
%draw_arc(sa,sl,su,colorAvgAng,offset,'-');
ax.RLim = [0 R];

% errors
sa = -39; % mean angle
sa=adjust-sa;
%if sa < 0, sa=sa+adjust2; end
su = 11.5;   % upper bound (smaller strike)
sl = 4.5;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
%ax=gca;
% draw_cone_hatched(sa-sl,sa+su,R,colorAvgAng);
ax.RLim = [0 R];

% M6, aftershocks of M6
colorAvgAng = 'b';
sa = 47.18; % mean angle
sa=adjust-sa;
su = 4;   % upper bound
sl = 3;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
% draw_line(sa,colorAvgAng,':',2);
offset=0.02;
%draw_arc(sa,sl,su,colorAvgAng,offset,':');
ax.RLim = [0 R];

% M6, strike from focal mechanism 
colorAvgAng = 'b';
sa = 49; % mean angle
sa=adjust-sa;
su = 2;   % upper bound
sl = 2;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
% draw_line(sa,colorAvgAng,'-',2);
offset=-0.01;
%draw_arc(sa,sl,su,colorAvgAng,offset,'-');
ax.RLim = [0 R];

% M6, strike from moment tensor 
colorAvgAng = 'b';
sa = 228-adjust2; % mean angle
su = 1;   % upper bound
sl = 1;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
% draw_line(sa,colorAvgAng,'--',2);
offset=-0.01;
%draw_arc(sa,sl,su,colorAvgAng,offset,'-');
ax.RLim = [0 R];

% errors
sa = 44.5; % mean angle
su = 5;   % upper bound
sl = 5;   % lower bound
%draw_cone(sa,sl,su,colorAvgAng);
% draw_cone_hatched(sa-sl,sa+su,R,colorAvgAng);


%text(deg2rad(88),2.3*F,'$\sigma_{Hmax}$','FontSize',14,'Rotation',88,...
%    'HorizontalAlignment','center',...
%    'VerticalAlignment','bottom', 'Interpreter','latex')
%  text(deg2rad(77),2.1*F,'$\dot{\varepsilon}_{Hmax}$','FontSize',16,'Rotation',86,...
%      'HorizontalAlignment','center','Color','m',...
%      'VerticalAlignment','bottom', 'Interpreter','latex')

% sigma 1 stress 
colorAvgAng = 'k';
sa =2.55;
sa=adjust-sa;
%draw_line(sa,colorAvgAng,'--',0.5);

% sigma 1 strain rate 
colorAvgAng = 'm';
sa =4.9;
sa=adjust-sa;
% draw_line(sa,colorAvgAng,'-',0.5);

ax.FontSize = 18;
ax.RTick = [];
ax.ThetaLim = [0 180];
%ax.ThetaLim = [270 90];
%ax.ThetaDir = 'clockwise';
%ax.ThetaZeroLocation = 'top';
ax.RLim = [0 R];

thetaticklabels({'East','60^\circ','30^\circ','North', '330^\circ','300^\circ','West'});

% Legend:
text(deg2rad(27),2.3*F,'left-lateral','FontSize',14,'Rotation',0,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','color','b')
%draw_h_line(37,1.5*F,'b','-','M6.4 focal mech.',2);
%draw_h_line(30,1.37*F,'b','--','M6.4 CMT',2);
%draw_h_line(21,1.28*F,'b',':','M6.4 aftershocks',2);
%draw_h_line(11,1.2*F,'b','none','active faults',2);

text(deg2rad(152),2.45*F,'right-lateral','FontSize',14,'Rotation',0,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','color','r')
%draw_h_line(164,2.83*F,'r','-','focal mech.',1);
%draw_h_line(169,2.78*F,'r','--','moment tensor',2);
%draw_h_line(174,2.73*F,'r',':','aftershocks',2);
%draw_h_line(160,2.9*F,'r','-','M7.1 focal mech.',2);
%draw_h_line(165,2.83*F,'r','--','M7.1 CMT',2);
%draw_h_line(170,2.79*F,'r',':','M7.1 aftershocks',2);
%draw_h_line(175,2.77*F,'r','','active faults');

%line ([1 2],[30 30],'lineStyle','-');
%saveas(gcf,'-r300','rose.png')
print(gcf,'-dpng','-r300','rose');
title("trans-tensional rose diagram")

function draw_line(sa,color,lineStyle,lineWidth)
  avgAngRad = deg2rad(sa);
  thetaAvgAng = [avgAngRad,avgAngRad];
  rL = rlim(gca);
  rhoAvgAng = rL;
  polarplot(gca,thetaAvgAng,rhoAvgAng,'lineStyle',lineStyle...
                  ,'lineWidth',lineWidth,'color',color), hold on

end

function draw_h_line(sa,r1,color,lineStyle,txt,wid)
  % draw a horizontal line on a polar plot starting at theta1,r1
  L=12; % length of a line
  theta1 = deg2rad(sa);
  theta2 = atan(r1*sin(theta1)/(r1*cos(theta1)+L));
  thetaAvgAng = [theta1,theta2];
  r2 = r1*sin(theta1)/sin(theta2);
  rhoAvgAng = [r1,r2];
  if length(lineStyle)>0
   polarplot(gca,thetaAvgAng,rhoAvgAng,'lineStyle',lineStyle...
            ,'lineWidth',wid,'color',color), hold on
  end

  L=2.1; % white space for a text label
  theta3 = atan(r2*sin(theta2)/(r2*cos(theta2)+L));
  r3 = r2*sin(theta2)/sin(theta3);
  text(theta3,r3,txt,'FontSize',12,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle')
end

function draw_cone(sa,sl,su,color)

  s1 = sa - sl;
  s2 = sa + su;
  sDist = deg2rad([s1 sa s2]);

  polarhistogram(sDist, 'BinLimits', [min(sDist) max(sDist)], 'FaceColor',color,'FaceAlpha',0.3), hold on
  avgAngRad = deg2rad(sa);
  thetaAvgAng = [avgAngRad,avgAngRad];
  rL = rlim(gca);
  rhoAvgAng = rL;
  polarplot(gca,thetaAvgAng,rhoAvgAng,'lineStyle','none'...
                  ,'lineWidth',2,'color',color), hold on
end

function draw_cone_hatched(s1,s2,R,color)

  L=2; % hatch interval 
  wid=0.5; % hatch line width
%  s1 = min(sDist(:));
%  s2 = max(sDist(:));
  theta1 = deg2rad(s1);
  theta2 = deg2rad(s2);
  trig=cos(theta1);
  rx1=R/trig;
  n=abs(floor(rx1/L));
  for i=1:n-1
   r1=abs(i*L/trig);
   r2=abs(i*L/cos(theta2));

   % draw a vert.  line on a polar plot starting at theta1,r1
   thetaAvgAng = [theta1,theta2];
   rhoAvgAng = [r1,r2];
   polarplot(gca,thetaAvgAng,rhoAvgAng,'lineStyle','-'...
                  ,'lineWidth',wid,'color',color), hold on
  end

%  rhoAvgAng = [0,2];
%size(rhoAvgAng)
%  polarplot(gca,thetaAvgAng,rhoAvgAng,'lineStyle','none'...
%                  ,'lineWidth',2,'color',color), hold on
end

function draw_arc(sa,sl,su,color,avgAngCiWhiskOffset,lineSp)

  s1 = sa - sl;
  s2 = sa + su;
  avgAngCi = sl+su;
  sDist = deg2rad([s1 sa s2]);
  avgAngRad = sDist(2);
  colorAvgAngCi = color;
  rL = rlim(gca);
  rhoAvgAng = rL;
  lineSp = '-';

  % distance between confidence-interval line and plot border in percent
  %  avgAngCiWhiskOffset = 0.02;
            avgAngCiLW = 2; % line width

            avgAngCiRad = deg2rad(avgAngCi);
            % ISNAN(AVGANGCI) == TRUE if the requirements for confidence levels are not
            % met, see CIRC_CONFMEAN line 70 -> do not plot

                self.avgAngCiH = gobjects(0);
                avgAngCiPlotArgs = {'lineStyle',lineSp,'color',colorAvgAngCi ...
                    ,'Clipping','off'};
                
                whiskWidthEnd = avgAngCiLW * 0.7; % width of whisker-endings
                thetaStepN = ceil(avgAngCi / 0.2); % plot line in 0.2-deg-steps
                if thetaStepN == 0, thetaStepN = 1; end
                thetaAvgAngCi = linspace(sDist(1),sDist(3),thetaStepN);
                rhoAvgAngCi = max(rL) + avgAngCiWhiskOffset * range(rL);
%                avgAngCiWhiskLen = abs(avgAngCiWhiskOffset)*2/3 * range(rL);
                avgAngCiWhiskLen = 0.02 * range(rL);
                rhoAvgAngCiWhiskEnd = ...
                    [rhoAvgAngCi + avgAngCiWhiskLen, rhoAvgAngCi - avgAngCiWhiskLen];
                rhoAvgAngCi = repmat(rhoAvgAngCi,thetaStepN,1);
                
                self.avgAngCiH(1,1) = polarplot(gca,thetaAvgAngCi,rhoAvgAngCi,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',avgAngCiLW,'Tag','avgAngCiWhisk');
                self.avgAngCiH(1,2) = polarplot(gca,[thetaAvgAngCi(1),thetaAvgAngCi(1)],rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                self.avgAngCiH(1,3) = polarplot(gca,[thetaAvgAngCi(end),thetaAvgAngCi(end)] ...
                    ,rhoAvgAngCiWhiskEnd,avgAngCiPlotArgs{:} ...
                    ,'lineWidth',whiskWidthEnd,'Tag','avgAngCiWhiskEnd');
                

end


