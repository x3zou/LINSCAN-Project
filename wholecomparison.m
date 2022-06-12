clear
%s: selected clusters
s1=importdata("selected_sqrt2.txt");
s2=importdata('selected_sqrt4.txt');
s3=importdata("selected_sqrt8.txt");
s4=importdata("selected_sqrt16.txt");

%l:clusters left after selection
l1=importdata("eps=sqrt4_input.txt");
l2=importdata("eps=sqrt8_input.txt");
l3=importdata("eps=sqrt16_input.txt");
l4=importdata("eps=sqrt16_remains.txt");

%entire catalog
real=importdata("realdata.txt");

figure();
line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s1(:,1),s1(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
line(l1(:,1),l1(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
title("eps=sqrt(2)")
legend("entire","chosen","remaining")

figure();
line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s2(:,1),s2(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
line(l2(:,1),l2(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
title("eps=sqrt(4)")
legend("entire","chosen","remaining")

figure();
line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s3(:,1),s3(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
line(l3(:,1),l3(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
title("eps=sqrt(8)")
legend("entire","chosen","remaining")

figure();
line(real(:,1),real(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s4(:,1),s4(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
line(l4(:,1),l4(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[1 0 1])
title("eps=sqrt(16)")
legend("entire","chosen","remaining")
