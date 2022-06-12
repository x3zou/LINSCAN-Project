clear
%a: all selected clusters
a1=importdata("correlation0.2iter_eps=sqrt2.txt");
a2=importdata("correlation0.2iter_eps=sqrt4.txt");
a3=importdata("correlation0.2iter_eps=sqrt8.txt");
a4=importdata("correlation0.2iter_eps=sqrt16.txt");
%s:clusters with further selection
s1=importdata("selected_sqrt2.txt");
s2=importdata('selected_sqrt4.txt');
s3=importdata("selected_sqrt8.txt");
s4=importdata("selected_sqrt16.txt");

figure();
line(a1(:,1),a1(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s1(:,1),s1(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
title("eps=sqrt(2)")
legend("linscan selected","further selected")

figure();
line(a2(:,1),a2(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s2(:,1),s2(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
title("eps=sqrt(4)")
legend("linscan selected","further selected")

figure();
line(a3(:,1),a3(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s3(:,1),s3(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
title("eps=sqrt(8)")
legend("linscan selected","further selected")

figure();
line(a4(:,1),a4(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5])
line(s4(:,1),s4(:,2),'lineStyle','none','Marker','.','MarkerSize',5,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1])
title("eps=sqrt(16)")
legend("linscan selected","further selected")