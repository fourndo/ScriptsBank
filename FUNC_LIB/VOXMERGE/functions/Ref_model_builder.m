% Build surface constraints using topocheck file and UBC mesh

close all
clear all

work_dir = 'C:\Projects\4080_Capstone_SaintRonans_grav_II\Modelling\Inversions\Input_25m\';

topo = load([work_dir 'topo_model.txt']);
geolfile =load([work_dir 'Geology_mean.den']);
Upbound_in =load([work_dir 'Geology_75pct.den']);
Lwbound_in =load([work_dir 'Geology_25pct.den']);

[mesh]=get_UBCmesh([work_dir 'Mesh_xy25m_z5m.msh']);

dx = mesh(3,mesh(3,:)~=0)';
dy = mesh(4,mesh(4,:)~=0)';
dz = mesh(5,mesh(5,:)~=0)';

nx = mesh(1,1); %size(X,1);    %number of cell in X
ny = mesh(1,2); %size(X,2);    %number of cell in Y
nz = mesh(1,3); %size(X,3);    %number of cell in Z

% Reshape 2D geological model and bounds
geolfile = reshape(geolfile,nx,ny);
Upbound_in = reshape(Upbound_in,nx,ny);
Lwbound_in = reshape(Lwbound_in,nx,ny);

nodatav = -100000;

topo = reshape(topo,nz,nx,ny);

% Allocate space for 3D models and assign no-data value
ref_model_surf = ones(nz,nx,ny) * nodatav;
ref_model_proj = ones(nz,nx,ny) * nodatav;
Upbound_out = ones(nz,nx,ny) * nodatav;
Lwbound_out = ones(nz,nx,ny) * nodatav;
Ws = ones(nz,nx,ny);
Wx = ones(nz,nx-1,ny)*5;
Wy = ones(nz,nx,ny-1)*5;
Wz = ones(nz-1,nx,ny)*5;

% Cycle through all the model columns and assign new ref-model and bounds
for jj = 1 : ny
    
    for ii = 1 : nx
        
        count = 0;
        
        for kk = 1 : nz
            
            
            
            if topo(kk,ii,jj) == 1 
                
                if count == 0
                    
                    top = kk;
                    
                end
                
                depth = sum(dz(top:kk))/dz(top);
                Ws(kk,ii,jj) = 10 /depth +1 ;
                
                if count < 2
%                     
%                    ref_model_surf(kk,ii,jj) = geolfile(ii,jj)-2.67;
%                    ref_model_proj(kk,ii,jj) = geolfile(ii,jj)-2.67;
%                    Upbound_out(kk,ii,jj) = Upbound_in(ii,jj)-2.67;
%                    Lwbound_out(kk,ii,jj) = Lwbound_in(ii,jj)-2.67;
%                    
%                    if Upbound_out(kk,ii,jj)==Lwbound_out(kk,ii,jj)
%                        
%                        Upbound_out(kk,ii,jj) = Upbound_out(kk,ii,jj)*1.1;
%                        Lwbound_out(kk,ii,jj) = Lwbound_out(kk,ii,jj) * 0.9;
%                        
%                    end
                   
%                     if ii<nx
% 
%                        Wx(kk,ii,jj) = 1;
% 
%                     end
%                    
%                     if jj<ny
% 
%                        Wy(kk,ii,jj) = 1;
% 
%                     end
                   
                   
                else
                    
%                    ref_model_surf(kk,ii,jj) = 0;

%                   if geolfile(ii,jj) > 2.98 && geolfile(ii,jj) < 3.0
%                       
%                       ref_model_proj(kk,ii,jj) = 2.75-2.67;
%                       
%                   else
%                       
%                       ref_model_proj(kk,ii,jj) = geolfile(ii,jj)-2.67;
%                       
%                   end
% 
%                 if geolfile(ii,jj) > 3.0 || kk<(-600)
%                    Lwbound_out(kk,ii,jj) = 0;
%                 else
%                    Lwbound_out(kk,ii,jj) = -0.5; 
%                 end
%                 
%                    Upbound_out(kk,ii,jj) = 1;
%                   
%                    
%                    if Upbound_out(kk,ii,jj)==Lwbound_out(kk,ii,jj)
%                        
%                        Upbound_out(kk,ii,jj) = Upbound_out(kk,ii,jj)*1.1;
%                        Lwbound_out(kk,ii,jj) = Lwbound_out(kk,ii,jj) * 0.9;
%                        
%                    end
                   
                end
                
                count = count+1;
                
            end
            
        end
    end
end

% ref_model_surf = reshape(ref_model_surf,nz*nx*ny,1);
ref_model_proj = reshape(ref_model_proj,nz*nx*ny,1);
Upbound_out = reshape(Upbound_out,nz*nx*ny,1);
Lwbound_out = reshape(Lwbound_out,nz*nx*ny,1);
Ws = reshape(Ws,nz*nx*ny,1);
Wx = reshape(Wx,nz*(nx-1)*ny,1);
Wy = reshape(Wy,nz*nx*(ny-1),1);
Wz = reshape(Wz,(nz-1)*nx*ny,1);

W = [Ws;Wx;Wy;Wz];
% Bounds = [Upbound_out Lwbound_out]
% save([work_dir 'ref_model_surf.dat'],'-ascii','ref_model_surf');
% save([work_dir 'ref_model_proj_v2.dat'],'-ascii','ref_model_proj');
% save([work_dir 'Upbound_out.dat'],'-ascii','Upbound_out');
% save([work_dir 'Lwbound_out.dat'],'-ascii','Lwbound_out');
save([work_dir 'Weights_v5.dat'],'-ascii','W');
       
                    
 