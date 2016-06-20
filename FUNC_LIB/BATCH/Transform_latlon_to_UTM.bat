cd C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\MtIsa\Data\DEM_AusGA

gdalwarp -s_srs EPSG:3112 -t_srs EPSG:32754 -of ERS DEM_AusGA.ers DEM_AusGA_UTM.ers
#gdalwarp -s_srs EPSG:4326 -t_srs EPSG:26910 -of ERS 092g06_0100_demw.dem 092g06_0100_demw_83UTM12.ers

# gdalwarp -s_srs EPSG:26911 -t_srs EPSG:26910 -of ERS 092h08_0100_deme_UTM11.ers 092h08_0100_deme_UTM10.ers
# gdalwarp -s_srs EPSG:26911 -t_srs EPSG:26910 -of ERS 092h08_0100_demw_UTM11.ers 092h08_0100_demw_UTM10.ers

pause