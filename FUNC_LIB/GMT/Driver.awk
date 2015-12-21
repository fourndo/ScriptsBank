BEGIN{


filin1="Obsfile.dat"
filin2="Predfile.dat"
ii=1

# Read down the observation file [X,Y,data]
while((getline<filin1)>0){
	X[ii]=$1
	Y[ii]=$2
	Obs[ii]=$3

	# Initialize [Xmin,Xmax,Ymin,Ymax]
	if(ii==1){
		Xmin=X[ii];
		Xmax=X[ii];
		Ymin=Y[ii];
		Ymax=Y[ii];
		Obsmin=Obs[ii];
		Obsmax=Obs[ii];
	}

	# Compare next entry to find max and min

	if(X[ii]<Xmin){
		Xmin=X[ii];
	}

	if(X[ii]>Xmax){
		Xmax=X[ii];
	}

	if(Y[ii]<Ymin){
		Ymin=Y[ii];
	}

	if(Y[ii]>Ymax){
		Ymax=Y[ii];
	}

	if(Obs[ii]>Obsmax){
		Obsmax=Obs[ii];
	}

	if(Obs[ii]>Obsmax){
		Obsmax=Obs[ii];
	}


	getline<filin2
	Pred=$3

	residual[ii]=Pred-Obs[ii]

	# Initialize residual [rmin,rmax]
	if(ii==1){
		rmin=residual[ii];
		rmax=residual[ii];
	}
	# Compare next entry to find [rmin,rmax]
	if(residual[ii]<rmin){
		rmin=residual[ii];
	}
	if(residual[ii]>rmax){
		rmax=residual[ii];
	}


	# Sum of residuals for std_dev
	sum+=residual[ii]

	# Print residual to file for GMT
	printf("% 9.4e\t% 9.4e\t% 8.5e\n",X[ii],Y[ii],residual[ii]) > "Residual.dat"
	ii++
	
}
 
nlines=ii-1

# Compute standard deviation 
for(ii=1;ii<=nlines;ii++){
	sumsq+=((residual[ii]-(sum/nlines))^2)
	std_dev=sqrt(sumsq/nlines)
}

# Print normalized residual
for(ii=1;ii<=nlines;ii++){
	rnorm[ii] = residual[ii]/std_dev

	# Initialize residual [rnormmin,rnormmax]
	if(ii==1){
		rnormmin=rnorm[ii];
		rnormmax=rnorm[ii];
	}
	# Compare next entry to find [rnormmin,rnormmax]
	if(rnorm[ii]<rnormmin){
		rnormmin=rnorm[ii];
	}
	if(rnorm[ii]>rnormmax){
		rnormmax=rnorm[ii];
	}

	printf("%9.4f\t%9.4f\t% 8.5e\n",X[ii],Y[ii],rnorm[ii]) > "Residual_norm.dat"
	
}


## Write observation color scales
# Round [min,max]
pow=0
while((Obsmax/(10^pow))>10){
	pow++
}

Obsmax= (int(Obsmax/(10^pow)) + 1) * (10^pow)

pow=0
while(sqrt((Obsmin/(10^pow))^2)>10){
	pow++
}

Obsmin= (int(Obsmin/(10^pow)) - 1) * (10^pow)

dr = int((Obsmax-Obsmin)/3)
line[1]=Obsmin" blue "(Obsmin+dr)" green"
line[2]=(Obsmin+dr)" green "(Obsmin+2*dr)" yellow"
line[3]=(Obsmin+2*dr)" yellow "(Obsmax)" red"

for(ii=1;ii<=3;ii++){
	printf("%s\n",line[ii]) > "colour.cpt"
}

delarray line

## Write residual color scales
# Round [min,max]
pow=0
while((rmax/(10^pow))>10){
	pow++
}

rmax= (int(rmax/(10^pow)) + 1) * (10^pow)

pow=0
while(sqrt((rmin/(10^pow))^2)>10){
	pow++
}

rmin= (int(rmin/(10^pow)) - 1) * (10^pow)

dr = int((rmax-rmin)/3)
line[1]=rmin" blue "(rmin+dr)" green"
line[2]=(rmin+dr)" green "(rmin+2*dr)" yellow"
line[3]=(rmin+2*dr)" yellow "(rmax)" red"

for(ii=1;ii<=3;ii++){
	printf("%s\n",line[ii]) > "colour_res.cpt"
}

delarray line

## Write normalized residual color scales
# Round [min,max]
#pow=0
#while((rnormmax/(10^pow))>10){
#	pow++
#}
#
#rnormmax= (int(rnormmax/(10^pow)) + 1) * (10^pow)
#
#pow=0
#while(sqrt((rnormmin/(10^pow))^2)>10){
#	pow++
#}
#
#rnormmin= (int(rnormmin/(10^pow)) - 1) * (10^pow)
rnormmax=3
rnormmin=-3
dr = int((rnormmax-rnormmin)/3)
line[1]=rnormmin" blue "(rnormmin+dr)" green"
line[2]=(rnormmin+dr)" green "(rnormmin+2*dr)" yellow"
line[3]=(rnormmin+2*dr)" yellow "(rnormmax)" red"

for(ii=1;ii<=3;ii++){
	printf("%s\n",line[ii]) > "colour_res_norm.cpt"
}

delarray line

# Compute ploting extents
trimx = 0.05 * (Xmax - Xmin)
trimy = 0.05 * (Ymax - Ymin)

pagexmin = int(Xmin - trimx)
pagexmax = int(Xmax + trimx)

xtick = int((pagexmax - pagexmin)/3)

pageymin = int(Ymin - trimy)
pageymax = int(Ymax + trimy)

ytick = int((pageymax - pageymin)/3)

## Mira Legend
line[1]="H 14 1 Mira Geoscience - AGIC"
line[2]="D 0 1p"
line[3]="L 10 0 C Title: Observed and Predicted data"
line[4]="D 0 1p"
line[5]="L 10 0 C Project: Synthetic example"
line[6]="D 0 1p"
line[7]="N 2"
line[8]="V 0 1p"
line[9]="L 10 0 C Scale: 1:10 000"
line[10]="L 10 0 C Date (YYYY-MM-DD): 2013-06-28"
line[11]="V 0 1p"
line[12]="D 0 1p"
line[13]="L 10 0 C Note: This is a test"
line[14]="G -0.15i"
line[15]="I Mira-Logo.eps 1i RT" 

for(ii=1;ii<=15;ii++){
	printf("%s\n",line[ii]) > "Mira.legend"
}


delarray line

## Write GMT script
line[1]="gmtset PAPER_MEDIA=11X17 LABEL_FONT_SIZE 10 HEADER_FONT_SIZE 16"
line[2]="gmtset PAGE_ORIENTATION=portrait ANNOT_FONT_SIZE_PRIMARY 10 HEADER_FONT_SIZE 12"
line[3]="gmtset BASEMAP_TYPE=plain HEADER_OFFSET 0.25c ANNOT_FONT_PRIMARY 0"
line[4]="gmtset COLOR_BACKGROUND=blue ANNOT_OFFSET_SECONDARY 0.05"
line[5]="gmtset MEASURE_UNIT=cm"
line[6]="del gridded_map.ps"
line[7]="psbasemap -R0/1.25/0/1.25 -JX20.5/29.25 -B:: -X0.25 -Y0.25 -K > gridded_map.ps"

line[8]="psscale -Ccolour.cpt -E -D6.1/18.0/6.0/0.25h -O -K >> gridded_map.ps"
line[9]="echo 0.35 0.79 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps"

line[10]="psscale -Ccolour.cpt -E -D14.8/18.0/6.0/0.25h -O -K >> gridded_map.ps"
line[11]="echo 0.86 0.79 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps"

line[12]="psscale -Ccolour_res.cpt -E -D14.8/5.0/6.0/0.25h -O -K >> gridded_map.ps"
line[13]="echo 0.86 0.235 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps"

line[14]="psscale -Ccolour_res_norm.cpt -E -D6.1/5.0/6.0/0.25h -O -K >> gridded_map.ps"
line[15]="echo 0.35 0.235 10 0 0 LT TMI(nT) | pstext -R -J -O -K -N >> gridded_map.ps"

line[16]="nearneighbor Obsfile.dat -R" pagexmin "/" pagexmax "/" pageymin "/" pageymax" -I100+ -S" 5*xtick " -GObsfile.grd "
line[17]="nearneighbor Predfile.dat -R" pagexmin "/" pagexmax "/" pageymin "/" pageymax" -I100+ -S" 5*xtick " -GPredfile.grd " 
line[18]="nearneighbor Residual.dat -R" pagexmin "/" pagexmax "/" pageymin "/" pageymax" -I100+ -S" 5*xtick " -GResidual.grd " 
line[19]="nearneighbor Residual_norm.dat -R" pagexmin "/" pagexmax "/" pageymin "/" pageymax" -I100+ -S" 5*xtick " -GResidual_norm.grd"

line[20]="grdgradient Obsfile.grd -GObs_surf_grad.grd -A0 -M -N"
line[21]="grdgradient Predfile.grd -GPred_surf_grad.grd -A0 -M -N"
line[22]="grdgradient Residual.grd -GResidual_surf_grad.grd -A0 -M -N"
line[23]="grdgradient Residual_norm.grd -GResidual_norm_surf_grad.grd -A0 -M -N"

# Plot Observed data
line[24]="grdimage Obsfile.grd -JX0.075m -IObs_surf_grad.grd -X2.6 -Y20.5 -Sn -Bg" ytick " -Bp" xtick ":\"\East (m)\":/" ytick ":.\"\Observed data\"::\"\North (m)\":WS --D_FORMAT=%%.0f  -Ccolour.cpt -O -K >> gridded_map.ps"
line[25]="grdcontour Obsfile.grd -R -J -A- -C500 -A20f14 -O -K -V >> gridded_map.ps"
#line[26]="psxy Obs_Model_Intrusive_2pc.dat -R -J -O -K -S+0.1 >> gridded_map.ps"
line[26]=""

# plot Predicted data
line[27]="grdimage Predfile.grd -JX0.075m -IPred_surf_grad.grd -X8.3 -Sn -Ccolour.cpt -O -K -Bg" ytick " -Bp" xtick ":\"\East (m)\":/" ytick ":.\"\Predicted data\":ES --D_FORMAT=%%.0f >> gridded_map.ps"
line[28]="grdcontour Predfile.grd -R -J -A- -C500 -A20f14 -O -K -V >> gridded_map.ps"
#line[29]="psxy maginv3d_010.pre -O -K -R -J -S+0.1 >> gridded_map.ps"
line[29]=""

# Plot Residual
line[30]="grdimage Residual.grd -JX0.075m -IResidual_surf_grad.grd -Y-13 -Sn -Ccolour_res.cpt -O -K -Bg" ytick " -Bp" xtick ":\"\East (m)\":/" ytick ":.\"\Residual\":ES --D_FORMAT=%%.0f >> gridded_map.ps"
line[31]="rem grdcontour Residual.grd -R -A- -J -C500 -A20f14 -O -K -V >> gridded_map.ps"
#line[32]="psxy Residual.dat -O -K -R -J -S+0.1  >> gridded_map.ps"
line[32]=""

# Plot Normalized Residual
line[33]="grdimage Residual_norm.grd -JX0.075m -IResidual_norm_surf_grad.grd  -X-8.3 -Sn -Ccolour_res_norm.cpt -O -K -Bg" ytick " -Bp" xtick ":\"\East (m)\":/" ytick ":.\"\Normalized Residual\"::\"\North (m)\":WS --D_FORMAT=%%.0f >> gridded_map.ps"
line[34]="rem grdcontour Residual_norm.grd -R -A- -J -C500 -A20f14 -O -K -V >> gridded_map.ps"
#line[35]="psxy Residual_norm.dat -O -K -R -J -S+0.1 >> gridded_map.ps"
line[35]=""

line[36]="pslegend -Dx7.65/-1.45i/8.075i/1.5i/TC -J -R -O -F Mira.legend -Gwhite >> gridded_map.ps"


for(ii=1;ii<=36;ii++){
printf("%s\n",line[ii]) > "GMT_plot.bat"
}




}