% ASCII to stats
% Load in ASCII files with various properties and rock units, and export
% statistics
% INPUT FILE FORMAT
% Prop1 Prop2 ... ID
% 
clear all 
close all

work_dir    = 'C:\LC\Private\dominiquef\Projects\4329_Goldcorp_Wabamisk_DCIP3D\Modelling\Inversion\3D';
prop_file   = 'Ref_model_v6.dat';
unit_file   = 'Geo_unit_model.dat';
outfile     = 'REF_STATS.dat';
unit_col = 6;
ndv = [1e-8 -1];

data = load([work_dir '\' prop_file]);
units = load([work_dir '\' unit_file]);

ID = unique(units);

fid = fopen([work_dir '\' outfile],'w');

for ii = 1 : length(ID)
    
    for jj = 1 : size(data,2)
        
        temp = data( units == ID(ii) & data(:,jj) ~= ndv(jj) , jj );
        lowb = prctile(temp,25);
        center = median(temp);
        uppb = prctile(temp,75);
        
        fprintf(fid,'< %3.1e ; %3.1e ; %3.1e >\t\t',lowb,center,uppb);
        
    end
    fprintf(fid,'%i',ID(ii));
    fprintf(fid,'\n');
    
end

fclose(fid);