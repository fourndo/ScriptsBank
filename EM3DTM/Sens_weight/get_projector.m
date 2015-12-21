function [P] = get_projector(Xc,Yc,Zc,iqx,iqy,iqz)
% Function P = get_projector(xc,yc,zc,qx,qy,qz)
% Create a sparse matrix to interpolate between querry cells
% The inportaltion uses a four-point inverse distance 
% weighting scheme. Assumes that Qx, Qy, Qz are evenly spaced querry
% locations using ndgrid
%
% Values for the entire mesh can be obtained by { y = P * x }, where 
% y is the full model and x are model values at the querry location
% [qx,qy,qz]
%
% INPUT
% xc,yc,zc : Coordinates of mesh cell center
% qx,qy,qz : Querry location inside mesh
nz = size(Xc,1); nx = size(Yc,2); ny = size(Zc,3);

nqx = length(iqx); nqy = length(iqy); nqz = length(iqz);

mcell = nz*nx*ny;
nq = nqz*nqx*nqy;

% Pre-allocate memory
P = spalloc(mcell,nq,4*nq);

fprintf('Start creating Projector matrix (P)\n');
fprintf('This process may take a while depending on the mesh size\n');

progress = -1;
tic
count = 1;

countj = 1;
for jj = 1 : ny
    
    if jj < iqy(1) || jj == iqy(countj) || jj > iqy(end)

        indy = countj;

    else

        indy = [countj countj+1];

        % Move to next querry node
        if jj > iqy(countj)
            countj = countj + 1;                
        end

    end
            
    counti = 1;    
    for ii = 1 : nx
                
            if ii < iqx(1) || ii == iqx(counti) || ii > iqx(end)
                
                indx = counti;
                
            else
                               
                indx = [counti counti+1];
                
                % Move to next querry node
                if ii > iqx(counti)
                    counti = counti + 1;                
                end
                
            end
        
        countk = 1;    
        for kk = 1 : nz
       
            if kk < iqz(1) || kk == iqz(countk) || kk > iqz(end)
                
                indz = countk;
                
            else
                               
                indz = [countk countk+1];
                
                % Move to next querry node
                if kk > iqz(countk)
                    countk = countk + 1;                
                end
                
            end
            
            [qqz,qqx,qqy] = ndgrid(indz,indx,indy);
            
            % Get index of selected querry cells
            tag = sub2ind([nqz nqx nqy],qqz(:),qqx(:),qqy(:));
            
            % Find closest cell
            nqcell = length(tag);

            P(count,tag) = 1/nqcell;

            d_iter = floor(count/mcell*20);

            count = count + 1;

            if  d_iter > progress

                fprintf('Processed %i pct of cells in %8.5f sec\n',d_iter*5,toc)
                progress = d_iter;

            end
            
        end
    end
end