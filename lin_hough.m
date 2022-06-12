clear all;
close all;

minlon = -60; maxlon = 25; minlat = 0; maxlat = 75; 
%minlon = -20; maxlon = -10; minlat = 60; maxlat = 72; 
%minlon = -10; maxlon = 10; minlat = -10; maxlat = 10; 
Fs=14;


fault_tr = 1; % flag for plotting the fault trace
% regional
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

%dx=0.125; dy=dx;
%xx=minlon:dx:maxlon; 
%yy=minlat:dy:maxlat; 
%[x,y]=meshgrid(xx,yy);
%z=0*x;

%for i=1:length(xx)
%  for j=1:length(yy)
%    for k=1:length(xlin(:,1))
%      if ((abs(xx(i)-xlin(k,1))<dx) & (abs(yy(j)-xlin(k,2))<dy))
%          z(j,i)=1.0;
%      end
%    end
%  end
%end
% save z z;

load z.mat;

figure(2)       
%pcolor(x,y,z), shading flat
pcolor(z), shading flat
c = gray;
c = flipud(c);
colormap(c);
colorbar
hold off

[H,theta,rho] = hough(z);

%figure(3)
imshow(imadjust(rescale(H)),[],...
       'XData',theta,...
       'YData',rho,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on
colormap(gca,hot)

  P = houghpeaks(H,25,'threshold',ceil(0.02*max(H(:))),'NHoodSize',[201 3] );

xp = theta(P(:,2));
yp = rho(P(:,1));
plot(xp,yp,'s','color','black'), hold off

lines = houghlines(z,theta,rho,P,'FillGap',20,'MinLength',2);

figure(4)
%imshow(z), hold on
  pcolor(z), shading flat, hold on
%return

max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
%plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
colormap(c);

