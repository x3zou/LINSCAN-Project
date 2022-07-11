
Instructions for running LINSCAN, selecting high quality clusters and making subplots

Xiaoyu Zou 07/05/2022 
x3zou@ucsd.edu

Before starting, please make sure you have all MATLAB and PYTHON subroutines

1.Download the earthquake catalog, load it into convert_ll_2_utm.m and run. Save the output xy coordinates and the converted catalog

2.Load the output xy coordinates into lin_scan_current (Andrew updated version) algorithm, save the output coordinates (temp), and temp_list which is used to identify each cluster. 

3.Go to HQselection.m, load the lin_scan input coordinates, the lin_scan output coordinates and the list file, run it, (set a suitable z value, initial value is 0, and please change it for iteration) and save the output ends, and clusters_saved file, which are the selected high-quality clusters, and their ends of best fitting lines.

(3.5.For iteration process, use iterationsteps.m to get the difference between the output high-quality clusters and the lin_scan input, and make it the new input for the lin_scan, repeat 2-3.5, then stack all the iterated data into one single file, which are the clusters_saved, and the ends, as input for next step. You may use outputstack.m)

4.Go to plotall.m, load the ends, the lin_scan input, the converted catalog, the high quality clusters (clusters_saved), the lin_scan output coordinates (maybe not necessary for now), the lin_scan list, and run it for subplots with clusters, best-fitting lines and composition focal mechanism. 

5.Go to lin_fm_new.m, load the same files as plotall.m, and plot the distribution of all focal mechanisms.Lin_fm_new also saves an output, "XXdihedral_input.txt",for the calculation of dihedral angle.

6. Go to dihedral_angle.m, load the "xxxdihedral_input.txt" from the previous step and run to get the dihedral angle histogram. Also go to rose_diag_mod.m, use load the same input to get the rose diagram.

(To make presentation plots, go to presentationplots.m)


