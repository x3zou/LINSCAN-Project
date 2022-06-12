%stack the iteration outputs
%Xiaoyu Zou, x3zou@ucsd.edu, 6/7/2022
%output:xx_ends.txt,xx_clusters_saved.txt
a=importdata("d2eps=2_ends.txt");
b=importdata("d2eps=4_ends.txt");
c=importdata("d2eps=8_ends.txt");
d=importdata("d2eps=16_ends.txt");
% e=importdata("d2eps=2_ends_2.txt");
f=importdata("d2eps=4_ends_2.txt");
g=importdata("d2eps=8_ends_2.txt");
% h=importdata("d2eps=16_ends_2.txt");
output=[a;b;c;d;f;g];
writematrix(output,"d2_ends.txt");
% 
% a=importdata("d2eps=2.txt");
% b=importdata("d2eps=4.txt");
% c=importdata("d2eps=8.txt");
% d=importdata("d2eps=16.txt");
% % e=importdata("d2eps=2_2.txt");
% f=importdata("d2eps=4_2.txt");
% g=importdata("d2eps=8_2.txt");
% % h=importdata("d2eps=16_2.txt");
% output=[a;b;c;d;f;g];
% writematrix(output,"d2_clusters_saved.txt");