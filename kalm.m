clear all;
close all;


x=[0:100]';
y=sin(x/10);

r=rand(size(x));
r=r-mean(r);
ind=find(x>70);
r(ind)=r(ind)*5;
z=y+r;

plot(x,z,'k*'), hold on
f=fit(x,z,'poly2')
plot(f,x,z)

return

%t=0:100;
%y=sin(y/10);
%t=x;
t = [0:100]';
u = sin(t/5);
n = length(t);
rng default
Q = 1;
R = 1;

w = sqrt(Q)*randn(n,1);
v = sqrt(R)*randn(n,1);
[out,x] = lsim(SimModel,[w,v,u]);

y = out(:,1);   % true response
ye = out(:,2);  % filtered response
yv = y + v;     % measured response
return
