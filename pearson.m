function [r] = pearson(nx,ny)
xa=mean(nx);
ya=mean(ny);
U=zeros(1,size(nx,1));
D1=zeros(1,size(nx,1));
D2=zeros(1,size(nx,1));
for k=1:size(nx,1)
    U(k)=(nx(k)-xa)*(ny(k)-ya);
    D1(k)=(nx(k)-xa)^2;
    D2(k)=(ny(k)-ya)^2;
end
r=sum(U)/sqrt(sum(D1)*sum(D2));
end
