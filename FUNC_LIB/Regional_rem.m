% Regional removal
clear all
close all

addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB';
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Documents\Projects\UBC_Tli_Kwi_Cho\Modelling\Inversion\Mag\Regional'

[H, I, Dazm, D, obsx, obsy, obsz, DIGHEM, wd] = read_MAG3D_obs([work_dir '\DIGHEM_Mag_2pc_floor9nt.obs']);
[H, I, Dazm, D, obsx, obsy, obsz, FWR_DIGHEM, wd] = read_MAG3D_obs([work_dir '\FWR_DIGHEM_regional.mag']);

DIGHEM_res = DIGHEM - FWR_DIGHEM;

write_MAG3D_TMI(work_dir,'DIGHEM_res_mag.obs',H,I,Dazm,...
    obsx,obsy,obsz,DIGHEM_res,wd)

[H, I, Dazm, D, obsx, obsy, obsz, Aerodat, wd_AD] = read_MAG3D_obs([work_dir '\aerotemData_noIGRF_uncert.mag']);
[H, I, Dazm, D, obsx, obsy, obsz, FWR_Aerodat, wd] = read_MAG3D_obs([work_dir '\FWR_Aerodat_regional.mag']);

Aerodat_res = Aerodat - FWR_Aerodat;

write_MAG3D_TMI(work_dir,'Aerodat_res_mag.obs',H,I,Dazm,...
    obsx,obsy,obsz,Aerodat_res,wd)
