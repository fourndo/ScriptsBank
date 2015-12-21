% Generate DCIP3D data location

clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\

work_dir = 'C:\LC\Private\dominiquef\Projects\4253_UBC_GIF_Testing\Codes\DCIP3D';

fid = fopen([work_dir '\DC_loc.loc'],'w');
% Write horizontal component
fprintf(fid,'!DC data\n');

xmin = -75;
xmax = 75;

xloc = xmin:10:xmax;
yloc = [-25 0 25];

for ii = 1 : length(yloc)
    
    for jj = 1 : length(xloc)
        
        fprintf(fid,'\n%12.5f %12.5f %12.5f %12.5f %i\n',xloc(jj),yloc(ii),xloc(jj),yloc(ii),length(xloc));
        for kk = 1 : length(xloc)
        fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',xloc(kk)-2.5,yloc(ii),xloc(kk)+2.5,yloc(ii),-99999,-99999);
        end
    end
    
end
        
 fclose(fid);       