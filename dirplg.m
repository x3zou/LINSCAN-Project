%subroutine for PTaxis function
function [dir,plg] = dirplg(e,n,u)
z=e*e+n*n;
z=sqrt(z);
plg=atan(-u/z);
if plg<0
    plg=-plg;
    e=-e;
    n=-n;
end
dir=atan(e/n);
end