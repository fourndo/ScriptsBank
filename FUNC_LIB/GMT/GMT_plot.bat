gmtset PAPER_MEDIA=11X17 LABEL_FONT_SIZE 10 HEADER_FONT_SIZE 16
gmtset PAGE_ORIENTATION=portrait ANNOT_FONT_SIZE_PRIMARY 10 HEADER_FONT_SIZE 12
gmtset BASEMAP_TYPE=plain HEADER_OFFSET 0.25c ANNOT_FONT_PRIMARY 0
gmtset COLOR_BACKGROUND=blue ANNOT_OFFSET_SECONDARY 0.05
gmtset MEASURE_UNIT=cm
del gridded_map.ps
psbasemap -R0/1.25/0/1.25 -JX20.5/29.25 -B:: -X0.25 -Y0.25 -K > gridded_map.ps
psscale -Ccolour.cpt -E -D6.1/18.0/6.0/0.25h -O -K >> gridded_map.ps
echo 0.35 0.79 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps
psscale -Ccolour.cpt -E -D14.8/18.0/6.0/0.25h -O -K >> gridded_map.ps
echo 0.86 0.79 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps
psscale -Ccolour_res.cpt -E -D14.8/5.0/6.0/0.25h -O -K >> gridded_map.ps
echo 0.86 0.235 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps
psscale -Ccolour_res_norm.cpt -E -D6.1/5.0/6.0/0.25h -O -K >> gridded_map.ps
echo 0.35 0.235 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps
nearneighbor Obsfile.dat -R385230/450695/5812197/5877662 -I100+ -S109105 -GObsfile.grd 
nearneighbor Predfile.dat -R385230/450695/5812197/5877662 -I100+ -S109105 -GPredfile.grd 
nearneighbor Residual.dat -R385230/450695/5812197/5877662 -I100+ -S109105 -GResidual.grd 
nearneighbor Residual_norm.dat -R385230/450695/5812197/5877662 -I100+ -S109105 -GResidual_norm.grd
grdgradient Obsfile.grd -GObs_surf_grad.grd -A0 -M -N
grdgradient Predfile.grd -GPred_surf_grad.grd -A0 -M -N
grdgradient Residual.grd -GResidual_surf_grad.grd -A0 -M -N
grdgradient Residual_norm.grd -GResidual_norm_surf_grad.grd -A0 -M -N
grdimage Obsfile.grd -JX0.075m -IObs_surf_grad.grd -X2.6 -Y20.5 -Sn -Bg21821 -Bp21821:"\East (m)":/21821:."\Observed data"::"\North (m)":WS --D_FORMAT=%%.0f  -Ccolour.cpt -O -K >> gridded_map.ps
grdcontour Obsfile.grd -R -J -A- -C500 -A20f14 -O -K -V >> gridded_map.ps

grdimage Predfile.grd -JX0.075m -IPred_surf_grad.grd -X8.3 -Sn -Ccolour.cpt -O -K -Bg21821 -Bp21821:"\East (m)":/21821:."\Predicted data":ES --D_FORMAT=%%.0f >> gridded_map.ps
grdcontour Predfile.grd -R -J -A- -C500 -A20f14 -O -K -V >> gridded_map.ps

grdimage Residual.grd -JX0.075m -IResidual_surf_grad.grd -Y-13 -Sn -Ccolour_res.cpt -O -K -Bg21821 -Bp21821:"\East (m)":/21821:."\Residual":ES --D_FORMAT=%%.0f >> gridded_map.ps
rem grdcontour Residual.grd -R -A- -J -C500 -A20f14 -O -K -V >> gridded_map.ps

grdimage Residual_norm.grd -JX0.075m -IResidual_norm_surf_grad.grd  -X-8.3 -Sn -Ccolour_res_norm.cpt -O -K -Bg21821 -Bp21821:"\East (m)":/21821:."\Normalized Residual"::"\North (m)":WS --D_FORMAT=%%.0f >> gridded_map.ps
rem grdcontour Residual_norm.grd -R -A- -J -C500 -A20f14 -O -K -V >> gridded_map.ps

pslegend -Dx7.65/-1.45i/8.075i/1.5i/TC -J -R -O -F Mira.legend -Gwhite >> gridded_map.ps
