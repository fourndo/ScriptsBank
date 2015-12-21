awk -f Obs_format.awk > Predfile.dat < maginv3d_004.pre 
awk -f Obs_format.awk > Obsfile.dat < 3986_Reg_200m.obs 

awk -f Driver.awk 

GMT_plot.bat