% WofE surface data

clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

%% INPUT FILES
work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\MtDore_files_for_Hypercube\Mt_Dore_Study_Hypercube_Expanded_Properties';
input_dir = [work_dir '\MtDore_2D_Targeting_Expanded_properties_UBC_format'];
out_dir = [work_dir '\Shell_properties'];

mesh = [work_dir '\MtDore_CEM_voxet.msh'];

% mfile{1} = 'MtDore_curvature_value_k_max_from_fault.dat';
% mfile{2} = 'MtDore_distance_to_high_fault_curvature_points.dat';
% mfile{3} = 'MtDore_subset_distance_to_faults.dat';
% mfile{4} = 'points_fault_curvature_k_max_at_intersection_with_faults.txt';
% mfile{5} = 'points_representing_interpolated_curvature_500m_from_faults.txt';
% mfile{2} = 'ERA_HighCurvature_Faults.mod';
% % mfile{7} = 'ERA_Dist_Rheological_Contact.mod';
% % mfile{8} = 'U_div_Th.mod';

% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\MtDore_files_for_Hypercube\Mt_Dore_Full3D_model\UBC_format_files\MtDore_gridded_grav_UBC_Format';
% mesh = [work_dir '\MtDore_Gridded_gravity.msh'];
% mfile{1} = 'Terrain_and_cover_corrected_gravity_data.mod';

% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\MtDore_files_for_Hypercube\Mt_Dore_Full3D_model\UBC_format_files\MtDore_gridded_TMI_UBC_Format';
% mesh = [work_dir '\MtDore_Gridded_TMI.msh'];
% mfile{1} = 'TMI.mod';
% 
% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\MtDore_files_for_Hypercube\Mt_Dore_Full3D_model\UBC_format_files\MtDore_Radiometrics_UBC_Format';
% mesh = [work_dir '\MtDore_Regional_Radiometrics.msh'];
% mfile{1} = 'Radiometrics_U_div_Th.mod';


%% Load the initial mesh and grab only the first cell
[xn,yn,zn] = read_UBC_mesh(mesh);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = ( xn(2:end) + xn(1:end-1) ) / 2;
yc = ( yn(2:end) + yn(1:end-1) ) / 2;
zc = ( zn(2:end) + zn(1:end-1) ) / 2;

mcell = nz * nx * ny;
[Zc,Xc,Yc] = ndgrid(zc,xc,yc);

% Look at first model and get nullcell
topocell = load([work_dir '\Topocells.dat']);
topocell = topocell==1;

cd(input_dir)
files_in = ls;

% Pre-allocate space for data matrix
datab = zeros(sum(topocell),36);

for ii = 1 : size(files_in,1)-2;
    
    temp = strtrim(files_in(ii+2,:));
    model = load([input_dir '\' temp]);

    % Reshape and propagate first cell to top
%     model = reshape(model,nz,nx,ny);
%     for kk = 1 : ny
%         
%         for jj = 1 : nx
%             
%             clmxy   = model(:,jj,kk);
%             grd     = nz - nnz(clmxy);
%             
%             if grd ~= nz
%                 clmxy(1:grd) = clmxy(grd+1)*ones(grd,1);
%             end
%             
%             model(:,jj,kk)=clmxy;
% 
% 
%         end
%         
%     end
%     
%     model = model(:);
%     model(model<0) = -99999;
    model(topocell & model==-99999) = 0;
    model(topocell == 0) = -99999;
    
    % save model out
    save([out_dir '\' temp(1:end-4) '_SHELL.dat'],'-ascii','model');
    
    prop{ii} = temp(1:end-4);
    datab(:,ii) = model(topocell);
end

fid = fopen([work_dir '\Mt_Dore_Database.dat'],'w');
% Write database to file
for ii  = 0 : size(datab,1)
    
    fprintf(fid,'%i\t',ii);
    
    for jj = 1 : size(datab,2)
        
        if ii == 0
            
            fprintf(fid,'%s\t',prop{jj});
            
        else
            
            fprintf(fid,'%12.8f\t',datab(ii,jj));
            
        end
        

    end
    
    fprintf(fid,'\n');
        
end

fclose(fid);
%% Load scatter points and interpolate on mesh
% kmax = load([work_dir '\' mfile{1}]);
% 
% index = kmax~=-99999 & topocell;
% kmax = kmax(index);
% kmax(kmax>0) = kmax(kmax>0) / max(kmax(kmax>0) );
% kmax(kmax<0) = kmax(kmax<0) / min(kmax(kmax<0) );
% % kmax = sqrt(kmax);
% % kmax(:,4) = abs(kmax(:,4));
% 
% % Pre-allocate space for model
% kmax_dist = zeros(nz,nx,ny);
% 
% for ii = 1 : nx
%         
%         for jj = 1 : ny
%                 
%             
%              R = (( abs( Xc(index) - xc(ii) ) + 25 ).^2 +...
%                         (abs(Yc(index) - yc(jj)) + 25 ).^2).^-0.5;
% 
% %             [r,index] = sort(R(:));
%             R = R + 1e-2 ;
% 
% %             for kk = 1 : 10
% %             wght = 1./R * kmax(index(kk),4);
% 
%             numer = sum( R .* kmax );
% 
%             denom = sum( R );
% 
% %             end
% 
%         kmax_dist(:,ii,jj) = numer./denom;
%    
%         end
% end
% 
% kmax_dist = kmax_dist(:);
% kmax_dist(topocell & kmax_dist==-99999) = 0;
% norm_kmax_dist = kmax_dist;
%     
% kmax_dist(topocell == 0) = -99999;
% save([work_dir '\' mfile{3}(1:end-4) '_Interp.dat'],'-ascii','kmax_dist');
% 
% % Normalize and dot multiply by distance to fault
% norm_kmax_dist(topocell) = kmax_dist(topocell)/ max(kmax_dist(topocell));
% 
% dist = load([work_dir '\' mfile{3}]);
% dist = 1./(dist).^0.25;
% dist = dist/min(dist);
% 
% norm_kmax_dist(topocell) = norm_kmax_dist(topocell) .* dist(topocell);
% norm_kmax_dist(topocell) = (norm_kmax_dist(topocell) ./  max(norm_kmax_dist(topocell))) ;
% norm_kmax_dist(topocell == 0) = -99999;
% save([work_dir '\lkmaxl_ldistl.dat'],'-ascii','norm_kmax_dist');

%% Load training site file and assign value no:0 or yes:1
% 
% tsites = load([work_dir '\..\MtDore_training_sites.txt']);
%    
% Xc = Xc(1,:,:);
% Yc = Yc(1,:,:);
% 
% % Training sites
% TS = zeros(nx,ny);
% 
% for ii = 1 : size(tsites,1)
%     
%     if tsites(ii,3) > yn(end) || tsites(ii,2) < xn(1) || tsites(ii,2) > xn(end) || tsites(ii,3) < yn(1) 
%         
%         continue
%         
%     end
%     
%     
%     % Find closest cell
%     R = sqrt( ( tsites(ii,2) - Xc ).^2 + ( tsites(ii,3) - Yc ).^2 );
%     
%     index = (R == min(min(R)));
%     
%     TS(index) = 1;
%     
% end
% 
% TS = TS==1;
% 
% topocell(:,TS) = topocell(:,TS) * 2;
% topocell = reshape(topocell,mcell,1);
% 
%  save([work_dir '\training.dat'],'-ascii','topocell');

