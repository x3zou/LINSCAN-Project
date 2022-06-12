%P and T axis calculation
function [azP, dipP, azT, dipT] = PTaxis(strike,dip,rake)
n=zeros(3);
s=zeros(3);
m=zeros(3);

ddir=strike+90;

z1=ddir*pi/180;
z2=dip*pi/180;
z3=rake*pi/180;

%normal vector to fault plane
n(1)=sin(z1)*sin(z2);
n(2)=cos(z1)*sin(z2);
n(3)=cos(z2);

%slickenside vector calculation
s(1)= -cos(z3)*cos(z1)-sin(z3)*sin(z1)*cos(z2);
s(2)= cos(z3)*sin(z1)-sin(z3)*cos(z1)*cos(z2);
s(3)= sin(z3)*sin(z2);

for i=1:2
    if i==1
        for j=1:3
            m(j)=s(j)-n(j);
        end
    [azP, dipP]=dirplg(m(1),m(2),m(3));
    else
        for j=1:3
            m(j)=s(j)+n(j);
        end
    [azT, dipT]=dirplg(m(1),m(2),m(3));
    end
end





