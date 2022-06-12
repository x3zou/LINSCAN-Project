function focalmech(fm, centerX, centerY, diam, varargin)
%	Draws full moment tensor beachball diagram of earthquake focal mechanism.
%      required inputs
%	fm: 3x3 or 1x6 vector (mrr, mtt, mff, mrt, mrf, mtf) or 
%          the six independent components of the moment tensor.
%	centerX, centerY: position to place beachball at position
%	diam: diameter of drawn beachball.
%
%	optional inputs
%	color for tensional region changed by 'r', 'g', 'b', 'y' or 
%       3 element color vector. default is black.
%   scalar argument is assumed to be an aspect ratio for non-equal axes.
%       Value is stretch in E-W direction.
%       'map' or 'Mercator' will do automatic Mercator aspect ratio
%       assuming centerY = lat
%   'dc' will force to best fitting double-couple solution
%   'nofill' will leave transparent (lines only)
%
%   Tensional and compressional regions are determined by sign of 
%   local 'radius'.
%   Negative length = compression, and positive length = tension
%   Length is found by linear combination of three eigen vectors.
%
%   Written to give MatLab equivalent of psmeca -Sm from GMT.
%   As such, this function requires a moment tensor. For planes (e.g.,
%   strikes and dips), use bb.m. bb is a nice beachball plotting function
%   from Oliver Boyd at CERI, but aimed strictly at double-couple
%   mechanisms.
%
%   Written by James Conder, Southern Illinois University, Jan 10, 2017
%   Published for use in
%   Conder, J.A. and C.A. Arciniegas, Fault reactivation in the Wabash  
%       Valley Fault Zone exhibited by the Bellmont 2012 earthquake, a Mt  
%       Carmel late aftershock, Seismological Research Letters,  
%       submitted, 2017

% defaults
fillcolor = [0 0 0];  % black for tensional region(s)
unitratio = 1;      % axis equal
DC = false;         % show non double-couple components
colorit = true;     % fill in tensional regions

% set parameters based on varargin inputs
if nargin > 4
    for i = 1:nargin-4
        if isnumeric(varargin{i})
            if isscalar(varargin{i})
                unitratio = varargin{i};
            end
            if length(varargin{i}) == 3
                fillcolor = varargin{i};
            end
        elseif ischar(varargin{i})
            switch varargin{i}
                case {'r'}
                    fillcolor = [ 1 0 0 ];
                case {'g'}
                    fillcolor = [ 0 1 0 ];
                case {'b'}
                    fillcolor = [ 0 0 1 ];
                case {'y'}
                    fillcolor = [ 1 1 0 ];
                case {'m'}
                    fillcolor = [ 1 0 1 ];
                case {'c'}
                    fillcolor = [ 0 1 1 ];
                
                case {'map','Map','mercator','Mercator'}
                    unitratio = 1/cos(centerY*pi/180);
                    
                case {'dc','DC','double-couple','doublecouple'}
                    DC = true;

                case {'nocolor','nofill','transparent'}
                    colorit = false;
            end 
        end
    end     
end

colormap([1 1 1; fillcolor])    % 2 color beachball fill 
%sE = unitratio*0.5*diam; sN = 0.5*diam; % scaling for plotting
sE = 0.5*diam; sN = 0.5*diam; % scaling for plotting
u = sE*cos(0:0.02:2*pi); w = sN*sin(0:0.02:2*pi);   % reference circle


%%% put fm (A&R convention, rtf) into 3x3 cartesian
%%% x = north, y = East, Z = down (Aki & Richards pg 113)
M = eye(3);
if length(fm) == 6
   M(1,1) = fm(2); M(2,2) = fm(3); M(3,3) = fm(1);
   M(2,1) = -fm(6); M(1,2) = M(2,1);
   M(3,1) = fm(4); M(1,3) = M(3,1);
   M(3,2) = -fm(5); M(2,3) = M(3,2);
else
   M(1,1) = fm(2,2); M(2,2) = fm(3,3); M(3,3) = fm(1,1);
   M(2,1) = -fm(3,2); M(1,2) = M(2,1);
   M(3,1) = fm(2,1); M(1,3) = M(3,1);
   M(3,2) = -fm(3,1); M(2,3) = M(3,2);
end
    
% find eigen values and eigenvectors of system and put in descending order
[V,D] = eig(M);
D = diag(D);
eig1 = max(D);      % largest eigenvalue
eig3 = min(D);      % smallest eigenvalue
if eig1 == eig3
    eig2 = eig1;
    i1 = 1; i2 = 2; i3 = 3;
else
    i1 = find(D == eig1,1);
    i3 = find(D == eig3,1);
    i2 = 6 - i1 - i3;
    eig2 = D(i2);
end

vT = V(:,i1);       % eigenvectors: TBP
vB = V(:,i2);
vP = V(:,i3);

% check for explosions or implosions
if eig1 < 0     % implosion
    fill(centerX+u,centerY+w,'w')   % fill white circle
    return
end
if eig3 > 0     % explosion
    fill(centerX+u,centerY+w,'k')   % fill black circle
    return
end

%%% if forced to be double couple
if DC
    eigB = 0.5*(eig1 - eig3);
    eig1 = eigB; eig2 = 0; eig3 = -eigB;
end

% Sweep out hemisphere to find negative & positive regions
dx = 0.02; dy = dx;     % grid loop
x = -1:dx:1; y = -1:dy:1;
uz = nan(length(x),length(y));  % lengths grid
vij = zeros(1,3);

for j = 1:length(y)
    yj = y(j);
    for i = 1:length(x)
        xi = x(i);
        trend = atan2(yj,xi);
        r2 = xi*xi + yj*yj;        
        plunge = pi/2 - 2*asin(sqrt(r2/2));  % equal area projection
        %plunge = pi/2 - 2*atan(sqrt(r2));  % equal angle projection
        
        if r2 <= 1       % inside beachball
            vij(1) = cos(trend)*cos(plunge);
            vij(2) = sin(trend)*cos(plunge);
            vij(3) = sin(plunge);
        
            % project vij onto principal stress directions to get weights
            uT = dot(vij,vT)*vT; wT = norm(uT)^2;
            uB = dot(vij,vB)*vB; wB = norm(uB)^2;  
            uP = dot(vij,vP)*vP; wP = norm(uP)^2;
        
            %%% find length (and sign) of vij
            uz(i,j) = wT*eig1 + wB*eig2 + wP*eig3;
        end
    end
end


%%% plot
% Note centerX & centerY are in conventional cartesian coords unlike x,y
% above where x is north.
hold on
if ~colorit
    plot(centerX+u,centerY+w,'k')   % plot circle boundary
    contour(centerX+sE*y,centerY+sN*x,uz,[0 0],'LineWidth',1);       % no fill
else
    fill(centerX+u,centerY+w,'w')   % fill white circle
    plot(centerX+u,centerY+w,'k')   % plot circle boundary


% fill in tensional regions
%contourf(centerX+sE*y,centerY+sN*x,uz,[0 0]);   % fill in tensional regions
% Need to use 'fill' as contourf will change pallette of entire figure

% get contours of zero level. Should be one or two segments
%c = contourc(y,x,uz,[0 0]);
[c,h] = contour(y,x,uz,[0 0]); delete(h);
c = c';
% Ugh. Would use contourc here instead of contour, but contourc matrix 
% puts uphill side on right instead of left (like contour & contourf) 


% find contour breaks
i = find(c(:,1) == 0);
j = find(c(:,2) >= 1);
i = intersect(i,j);
if ~isempty(find(c(i,2) <= 3, 1)) % remove contours with 3 or fewer points
    for j = fliplr(i')
        if c(j,2) <=3
            c(j:j+c(j,2),:) = [];
        end
    end
    i = find(c(:,2) >= 1);
end

% close and fill
if length(i) == 1       % 1 closed contour
    c(1,:) = [];        % remove contour header 
    % closed contour must have vertical plunge in it (check for sign)
    vij = [ 0 0 1];
    uT = dot(vij,vT)*vT; wT = norm(uT)^2;
    uB = dot(vij,vB)*vB; wB = norm(uB)^2;  
    uP = dot(vij,vP)*vP; wP = norm(uP)^2;
        
    %%% find length (and sign) of vij
    uv = wT*eig1 + wB*eig2 + wP*eig3;
    
    if uv > 0
        fill(centerX+sE*c(:,1),centerY+sN*c(:,2),fillcolor)
    else 
        fill(centerX+u,centerY+w,fillcolor)
        fill(centerX+sE*c(:,1),centerY+sN*c(:,2),'w')
    end
else                    % 2 unclosed contours
    % find 4 spots where contours touch edge
    ia = ([ 2 i(2)-1 i(2)+1 length(c)])';
    ang = atan2(c(ia,2),c(ia,1));
    
    % force angle to increase ccw from angle2
    dshift = ang(2);
    ang = ang - dshift;
    j = find(ang < 0);
    ang(j) = ang(j) + 2*pi;
    ang = ang + dshift;
    
    if (ang(1) - ang(2)) < pi      % close contours on themselves
        a1 = (ang(2):0.02:ang(1)-0.02)';    % arc from 2 to 1
        ac = [cos(a1) sin(a1)];
        c1 = [c(2:i(2)-1,:); ac; c(2,:)];
        fill(centerX+sE*c1(:,1),centerY+sN*c1(:,2),fillcolor)
        
        a1 = (ang(4):0.02:ang(3)-0.02)';
        ac = [cos(a1) sin(a1)];
        c1 = [c(i(2)+1:end,:); ac; c(i(2)+1,:)];
        fill(centerX+sE*c1(:,1),centerY+sN*c1(:,2),fillcolor)
    else
        % connect end of 1 to start of 2 along boundary and
        % end of 2 to start of 1

        a1 = (ang(2):0.02:ang(3)-0.02)';
        ac1 = [cos(a1) sin(a1)];
        a2 = (ang(4):0.02:ang(1)-0.02)';
        ac2 = [cos(a2) sin(a2)];
        c1 = [c(2:i(2)-1,:); ac1; c(i(2)+1:end,:); ac2; c(2,:)];
        fill(centerX+sE*c1(:,1),centerY+sN*c1(:,2),fillcolor)
    end
end
end
        
    
end

        
