cd C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\GIF_model\092g06

gdalwarp -s_srs EPSG:4326 -t_srs EPSG:26910 -of ERS 092g06_0100_deme.dem 092g06_0100_deme_83UTM12.ers
gdalwarp -s_srs EPSG:4326 -t_srs EPSG:26910 -of ERS 092g06_0100_demw.dem 092g06_0100_demw_83UTM12.ers

# gdalwarp -s_srs EPSG:26911 -t_srs EPSG:26910 -of ERS 092h08_0100_deme_UTM11.ers 092h08_0100_deme_UTM10.ers
# gdalwarp -s_srs EPSG:26911 -t_srs EPSG:26910 -of ERS 092h08_0100_demw_UTM11.ers 092h08_0100_demw_UTM10.ers

pause