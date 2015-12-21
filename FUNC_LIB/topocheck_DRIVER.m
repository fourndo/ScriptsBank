%% Topocheck_DRIVER
clear all
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\HyperCube\Processing\MtDore_files_for_Hypercube\Targeting_2015';

meshfile = 'Mesh_100m.msh';
topofile = 'Topo_SRTM_lim.topo';
dist_w = 'dist_weight.txt';

ndv = 1e-20;

%% SCRIPT START HERE
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

nx = ( length(xn) - 1 );
ny = ( length(yn) - 1 );
nz = ( length(zn) - 1 );

m = load([work_dir '\' dist_w]);

m = reshape(m,nz,nx,ny);
m(m~=ndv) = 1;
m(m==ndv) = 0;

for ii = 1 : ny
    
    for jj = 1 : nx
        
        for kk = 2 : nz
            
            if m(kk-1,jj,ii)==0 && m(kk,jj,ii)==1 
                
                m(kk:end,jj,ii) = 0;
                m(kk,jj,ii) = 1;
                
                break
                
            end
            
        end

    end
    
end

act_lay = m(:);

% [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
% 
% mcell = nx * ny * nz;
% 
% if isempty(topofile)==1
%     
%     nxn = length(xn);
%     nyn = length(yn);
%     topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
%     
% else
%     % Load topo
%     topo = read_UBC_topo([work_dir '\' topofile]);
% end
% 
% % Create nullcell
% [nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% 
% % Output active layer
% act_lay = nullcell*-1;
% 
% act_lay(tcellID) = 1;

save([work_dir '\active_cell.dat'],'-ascii','act_lay');
