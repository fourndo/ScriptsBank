cd C:\LC\Private\dominiquef\Projects\LivingstonCreek_VTEM\Data\DEM

gdalwarp -s_srs EPSG:4617 -t_srs EPSG:32608 -of ERS 105e01_0100_deme.dem 105e01_0100_deme.ers
gdalwarp -s_srs EPSG:4617 -t_srs EPSG:32608 -of ERS 105e01_0100_demw.dem 105e01_0100_demw.ers

pause