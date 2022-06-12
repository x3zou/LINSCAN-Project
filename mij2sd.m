function [S1,D1,S2,D2]=mij2sd(fm)
%function [S1,D1,S2,D2]=mij2sd(fm)
%
%	caldulates strike and dip S1, D1, S2, D2
%		of the focal planes. Can be vectors of
%		multiple focal mechanisms.
%	fm - focal mechanism that is  number of mechnisms (NM) 
%		NM x 6 (mxx, myy, mzz, mxy, mxz, myz - the six independent components of
%		the moment tensor). The strike is of the first plane, clockwise relative to north.
%		The dip is of the first plane, defined clockwise and perpedicular to strike, 
%		relative to horizontal such that 0 is horizontal and 90 is vertical. The rake is 
%		of the first focal plane solution. 90 moves the hanging wall up-dip (thrust),
%		0 moves it in the strike direction (left-lateral), -90 moves it down-dip
%		(normal), and 180 moves it opposite to strike (right-lateral).
%


[ne,n] = size(fm);
if n == 6
	for j = 1:ne
		[s1(j),d1(j),r1(j)] = mij2sdr(fm(j,1),fm(j,2),fm(j,3),fm(j,4),fm(j,5),fm(j,6));
	end
else
	s1 = fm(:,1);
	d1 = fm(:,2);
	r1 = fm(:,3);
end

r2d = 180/pi;
d2r = pi/180;


mech = zeros(ne,1);
j = find(r1 > 180);
r1(j) = r1(j) - 180;
mech(j) = 1;
j = find(r1 < 0);
r1(j) = r1(j) + 180;
mech(j) = 1;

% Get azimuth and dip of second plane
[s2,d2,r2] = AuxPlane(s1,d1,r1);

for ev = 1:ne
	S1 = s1(ev);
	D1 = d1(ev);
	S2 = s2(ev);
	D2 = d2(ev);
	P = r1(ev);

if D1 >= 90
   D1 = 89.9999;
end
if D2 >= 90
   D2 = 89.9999;
end

end

function [strike, dip, rake] = AuxPlane(s1,d1,r1);
%function [strike, dip, rake] = AuxPlane(s1,d1,r1);
% Get Strike and dip of second plane, adapted from Andy Michael bothplanes.c
r2d = 180/pi;

z = (s1+90)/r2d;
z2 = d1/r2d;
z3 = r1/r2d;
%/* slick vector in plane 1 */
sl1 = -cos(z3).*cos(z)-sin(z3).*sin(z).*cos(z2);
sl2 = cos(z3).*sin(z)-sin(z3).*cos(z).*cos(z2);
sl3 = sin(z3).*sin(z2);
[strike, dip] = strikedip(sl2,sl1,sl3);

n1 = sin(z).*sin(z2);  %/* normal vector to plane 1 */
n2 = cos(z).*sin(z2);
n3 = cos(z2);
h1 = -sl2; %/* strike vector of plane 2 */
h2 = sl1;
%/* note h3=0 always so we leave it out */

z = h1.*n1 + h2.*n2;
z = z./sqrt(h1.*h1 + h2.*h2);
z = acos(z);

rake = zeros(size(strike));
j = find(sl3 > 0);
rake(j) = z(j)*r2d;
j = find(sl3 <= 0);
rake(j) = -z(j)*r2d;

function [strike, dip] = strikedip(n, e, u)
%function [strike, dip] = strikedip(n, e, u)
%       Finds strike and dip of plane given normal vector having components n, e, and u
%

% Adapted from Andy Michaels stridip.c

r2d = 180/pi;

j = find(u < 0);
n(j) = -n(j);
e(j) = -e(j);
u(j) = -u(j);

strike = atan2(e,n)*r2d;
strike = strike - 90;
while strike >= 360
        strike = strike - 360;
end
while strike < 0
        strike = strike + 360;
end

x = sqrt(n.^2 + e.^2);
dip = atan2(x,u)*r2d;


function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
%function [str,dip,rake] = mij2sdr(mxx,myy,mzz,mxy,mxz,myz)
%
%   INPUT
%	mij - siz independent components of the moment tensor
%
%   OUTPUT
%	str - strike of first focal plane (degrees)
%	dip - dip of first focal plane (degrees)
%	rake - rake of first focal plane (degrees)
%
%Adapted from code, mij2d.f, created by Chen Ji and given to me by Gaven Hayes.
%

a = [mxx mxy mxz; mxy myy myz; mxz myz mzz];
[V,d] = eig(a);

D = [d(3,3) d(1,1) d(2,2)];
V(2:3,1:3) = -V(2:3,1:3);
V = [V(2,3) V(2,1) V(2,2); V(3,3) V(3,1) V(3,2); V(1,3) V(1,1) V(1,2)];

IMAX = find(D == max(D));
IMIN = find(D == min(D));
AE = (V(:,IMAX)+V(:,IMIN))/sqrt(2.0);
AN = (V(:,IMAX)-V(:,IMIN))/sqrt(2.0);
AER = sqrt(AE(1)^2+AE(2)^2+AE(3)^2);
ANR = sqrt(AN(1)^2+AN(2)^2+AN(3)^2);
AE = AE/AER;
AN = AN/ANR;
if (AN(3) <= 0.)
	AN1 = AN;
	AE1 = AE;
else
	AN1 = -AN;
	AE1 = -AE;
end 
[ft,fd,fl] = TDL(AN1,AE1);
str = 360 - ft;
dip = fd;
rake = 180 - fl;

function [FT,FD,FL] = TDL(AN,BN)
XN=AN(1);
YN=AN(2);
ZN=AN(3);
XE=BN(1);
YE=BN(2);
ZE=BN(3);
AAA=1.0E-06;
CON=57.2957795;
if (abs(ZN) < AAA)
	FD=90.;
	AXN=abs(XN);
	if (AXN > 1.0) 
		AXN=1.0;
	end
	FT=asin(AXN)*CON;
	ST=-XN;
	CT=YN;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	FL=asin(abs(ZE))*CON;
	SL=-ZE;
	if (abs(XN) < AAA) 
		CL=XE/YN;
	else
		CL=-YE/XN;
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
else
	if (-ZN > 1.0) 
		ZN=-1.0;
	end
	FDH=acos(-ZN);
	FD=FDH*CON;
	SD=sin(FDH);
	if  (SD == 0)
		return;
	end 
	ST=-XN/SD;
	CT=YN/SD;
	SX=abs(ST);
	if (SX > 1.0) 
		SX=1.0;
	end
	FT=asin(SX)*CON;
	if (ST >= 0. & CT < 0) 
		FT=180.-FT;
	end
	if (ST < 0. & CT <= 0) 
		FT=180.+FT;
	end
	if (ST < 0. & CT > 0) 
		FT=360.-FT;
	end
	SL=-ZE/SD;
	SX=abs(SL);
	if (SX > 1.0) 
		SX=1.0;
	end
	FL=asin(SX)*CON;
	if (ST == 0) 
		CL=XE/CT;
	else
		XXX=YN*ZN*ZE/SD/SD+YE;
		CL=-SD*XXX/XN;
		if (CT == 0) 
			CL=YE/ST;
		end
	end 
	if (SL >= 0. & CL < 0) 
		FL=180.-FL;
	end
	if (SL < 0. & CL <= 0) 
		FL=FL-180.;
	end
	if (SL < 0. & CL > 0) 
		FL=-FL;
	end
end 
