cd C:\Projects\4180_Wallbridge_Wisner_IP\Export
ogr2ogr -f DXF -overwrite -s_srs EPSG:26717 -t_srs EPSG:4326 4180_CDED_countour_20m.dxf 
ogr2ogr -f "KML" 4180_CDED_countour_20m_LONLAT.kml 4180_CDED_countour_20m_LONLAT.dxf
pause