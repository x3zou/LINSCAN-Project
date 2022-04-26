function [kld] = getkld(data)
[x1,y1,x2,y2]=lin_fit2(data(:,1),data(:,2));
vector=[x1,y1;x2,y2];
projpoint=zeros(1,2);
newpoint=zeros(1,size(data,1));
for i=1:size(data,1)
    pj=proj(vector,data(i,:));
    pj=pj';
    projpoint=vertcat(projpoint,pj);
end
projpoint(1,:)=[];
% line([x1 x2],[y1 y2],'lineStyle','-','Color','b','LineWidth',1);
% hold on
% scatter(projpoint(:,1),projpoint(:,2))
x0=min(x1,x2);
if x0==x1
    y0=y1;
else
    y0=y2;
end
for i=1:size(projpoint,1)
    if projpoint(i,1)>x0
        newpoint(1,i)=sqrt((projpoint(i,1)-x0)^2+(projpoint(i,2)-y0)^2);
    else
        newpoint(1,i)=-sqrt((projpoint(i,1)-x0)^2+(projpoint(i,2)-y0)^2);
    end
end
% pd1=fitdist(newpoint','kernel');
h=histogram(newpoint,'Normalization','probability');
set(h,'Visible','off')
what=h.Values;
% ywhat=pdf(pd1,newpoint);
u=ones(1,size(what,2))*1/size(what,2);
kld=getKullbackLeibler(u,what);
end