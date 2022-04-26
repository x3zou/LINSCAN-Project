real=importdata('realdata.txt');
ends=importdata('endtotal.txt');
line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',8,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
hold on
for i=1:size(ends,1)
    x1=ends(i,1);
    y1=ends(i,2);
    x2=ends(i,3);
    y2=ends(i,4);
    line([x1 x2],[y1 y2],'lineStyle','-','Color','b','LineWidth',1)
end
