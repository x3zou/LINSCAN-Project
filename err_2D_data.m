function [Err] = err_2D_data(k, XData, YData, P)
%
% function [Err] = err_2D_data(XData, YData, P)
%
% Orthogonal linear regression method in 2D for model: y = a + bx   
%
% Input parameters:
%  - XData: input data block -- x: axis
%  - YData: input data block -- y: axis
%  - P: optional vector of model parameters [b-slope, a-offset]
%       if P is empty, compute it using robust least squares
%
% Return parameters:
%  - Err: error in the best-fit slope estimate (b) 
%
%
% YF Date: 08/14/2020
%
kx=length(XData);
ky=length(YData);
if kx ~= ky
   disp('Incompatible X and Y data.');
   close all;
end

%%if length(P(:) >0)
% Yhat = XData.*P(1) + P(2);
% Xhat = ((YData-P(2))./P(1));

% alpha = atan(abs((Yhat-YData)./(Xhat-XData)));
%%else
%[p1,s1]=polyfit(XData,YData,1);
%[py,sy]=polyfit(YData,XData,1);

%s0=sqrt(sum((Yhat(:)-mean(Yhat(:))).^2)/sum((Xhat(:)-mean(Yhat(:))).^2));
%s0=sqrt(sum((YData(:)-mean(YData(:))).^2)/sum((XData(:)-mean(YData(:))).^2));

s0=P(1);


%yval=p1(2)+p1(1)*XData;
%[yval,delta]=polyval(p1,XData,s1);

mdlr = fitlm(XData,YData,'RobustOpts','on');
beta1 = mdlr.Coefficients.Estimate;
p1(2)=beta1(1); p1(1)=beta1(2); 

yval=XData*p1(1)+p1(2);
%line(XData,yval,'LineStyle','-','Color','m','Linewidth',1.5), hold on


%[yval,delta]=polyval(p1,XData,s1);
mdlr = fitlm(YData,XData,'RobustOpts','on');
beta2 = mdlr.Coefficients.Estimate;
py(2)=beta2(1); py(1)=beta2(2); 

p2(1)=1/py(1);
p2(2)=-py(2)/py(1);

yval=XData*p2(1)+p2(2);
%line(XData,yval,'LineStyle','-','Color','b','Linewidth',1.5), hold on

% slope estimates:
P1=p1(1);  P2=p2(1);
fact=1.18;
fact=3;

%if k==11, fprintf('%f %f \n',P1,P2), end

difang = abs(atan(P1)-atan(P2))/pi*180;
if abs(atan(P1)-atan(P2)) > pi/3 % no good fit found
  P1=s0+tan(abs(atan(P1))/fact);
  P2=s0+tan(abs(atan(P2))/fact);
  fprintf('no acceptable fit found for cluster %d  \n',k);
end
%if k==11, fprintf('%f %f \n',P1,P2), end

thresh=20; % threshold for max. slope; rely on the other estimate
if abs(P1) > thresh
 P1=P2;
end
if abs(P2) > thresh
 P2=P1;
end

if abs(atan(P1)-atan(s0)) > pi/3
  P1=s0+tan(abs(atan(P1))/fact);
  fprintf('slope P1 constrained for cluster %d  \n',k);
end

if abs(atan(P2)-atan(s0)) > pi/3
  P2=s0+tan(abs(atan(P2))/fact);
  fprintf('slope P2 constrained for cluster %d  \n',k);
end
%if k==11, fprintf('%f %f \n',P1,P2), end

%if max([abs(atan(p1(1))-atan(s0)) abs(atan(p2(1))-atan(s0))]) > pi/4  
% P1=s0+tan(abs(atan(p1(1)) - atan(p2(1))));
% P2=P1;
% k
%end

%if abs(atan(p2(1))-atan(s0)) > pi/4
%%  P2 =p1(1);
%P2 =tan(mean(atan([p1(1) p2(1)])));
%k
%else
%  P2=p2(1);
%end

%d0=(abs(abs(p2(1))-s0)+abs(abs(p1(1))-s0))/2; % strike error, in rad
%d0=(abs(abs(P2)-s0)+abs(abs(P1)-s0))/2; % strike error, in rad

d0=(abs(P2-s0)+abs(P1-s0))/2; % strike error, in rad

%  fprintf('%d %f %f %f  %f \n', k, s0, p1(1), p2(1), d0); 
%  fprintf('%d %f %f %f  %f \n', k, s0, P1, P2, d0); 

d1=atan(s0-d0)/pi*180;
d2=atan(s0+d0)/pi*180;
%Err=abs(d2-d1)/2; % strike error, in deg. 
Err=abs(d2-d1); % strike error, in deg. 

fact1=4;
if Err < difang/fact1
 Err=difang/fact1;
 fprintf('Error estimate increased for cluster %d  \n',k);
end

%fprintf('%d %f %f %f %f  %f \n', k, atan(s0)/pi*180, d1, d2, difang/fact1, Err); 

%d=abs(Xhat-XData).*sin(alpha);
%%Err=sum(abs(d));
%Err=sum(d.^2);
