function [OMESid] = MAG3C_get_OMES_id(nullcell,nlayer)
% Function MAG3C_get_OMES_id(nullcell,acellID,1)
% Return the cell ID of the top most active cells for nlayers
%
% INPUT: 
% nullcell: Binary matrix of 1-0 for active-inactive cells
%           DIM(nz x nx x ny)
% nlayer: Integer for number of layers to be used below topography
%
% OUTPUT
% OMESid: cell id number of active cells within the OMES layer

nz = size(nullcell,1);
nx = size(nullcell,2);
ny = size(nullcell,3);

OMESid = zeros(nx*ny*nlayer,1);

countID = 1;
for jj = 1 : ny
    
    for ii = 1 : nx
        
        countn = 1;
        for kk = 1 : nz
        
            if nullcell(kk,ii,jj) == 0
                
                continue
                
            else
               
                OMESid(countID) = sub2ind([nz nx ny],kk,ii,jj);
                countID = countID + 1;
                countn = countn + 1;
            end
            
            if countn > nlayer
                
                break
                
            end
            
        end
                
    end
    
end