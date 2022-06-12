function [fm]=sdr2mij(v3,deg2rad)
%function sdr2mij(v3,deg2rad)
%
%	caldulates moment tensor given v3 = strike, dip, and rake
%       Can be vectors of
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
%       deg2rad - 1 if already in radians, pi/180 if in degrees
%convert to radians:
s=v3(1)*deg2rad;
d=v3(2)*deg2rad;
r=v3(3)*deg2rad;

fm(2)=-(sin(d)*cos(r)*sin(2*s)+sin(2*d)*sin(r)*(sin(s))^2);
fm(3)=(sin(d)*cos(r)*sin(2*s)-sin(2*d)*sin(r)*(cos(s)^2));
fm(1)=sin(2*d)*sin(r);
fm(6)=-(sin(d)*cos(r)*cos(2*s)+0.5*sin(2*d)*sin(r)*sin(2*s));
fm(4)=-(cos(d)*cos(r)*cos(s)+cos(2*d)*sin(r)*sin(s));
fm(5)=(cos(d)*cos(r)*sin(s)-cos(2*d)*sin(r)*cos(s));

end 
