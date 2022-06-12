%Get the input for next step of iteration and save all necessary variables
%please run this immediately after HQselection.m, do not clear variables
%!!!please remember to update the file names!!!
%Xiaoyu Zou, x3zou@ucsd.edu, 5/31/2022
writematrix(clusters_saved,"d1eps=16_2.txt")
writematrix(ends,"d2eps=8_ends_2.txt")
diff=real(~ismember(real,points_added,'rows'),:);
writematrix(diff,"d1eps=16_output_2.txt")%input for the next round of iteration

