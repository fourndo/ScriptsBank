function [Tx,Ty,Tz] = MAG3C_T_octree( Tx,Ty,Tz,Xnode,Ynode,Znode,nullcell,topocell,toponode,sprcell,obsx,obsy,obsz)
% [Tx,Ty,Tz] = MAG3C_T_octree( Tx,Ty,Tz, nullcell,topocell,toponode,dobsx,dobsy,dobsz)

% [idz,idx,idy] = ndgrid(1:length(dobsz)-1,1:length(dobsx)-1,1:length(dobsy)-1);
nz = size(toponode,1)-1;
nx = size(toponode,2)-1;
ny = size(toponode,3)-1;

% Need index within 
mcell = nx*ny*nz;
tij = [Tx;Ty;Tz];
tsuper = tij;
for jj = 1 : length(topocell)
    
    cellID = topocell(jj);
    superC = sprcell(jj);
    
    [K,I,J] = ind2sub([nz,nx,ny],cellID);
            
    ttemp = tsuper;
    
    dx = Xnode(K,I+1,J) - Xnode(K,I,J);
    dy = Ynode(K,I,J+1) - Ynode(K,I,J);
    dz = Znode(K,I,J) - Znode(K+1,I,J);
    
%     ztop = toponode(1,I:I+1,J:J+1);
    
    % Reshape topography for the current cell (4 points)
    topo = [reshape(Xnode(1,I:I+1,J:J+1),4,1) reshape(Ynode(1,I:I+1,J:J+1),4,1) reshape(toponode(1,I:I+1,J:J+1),4,1)];
        
    % Breakdown the topocell until kernel does not change more than 1%
    dtij = 1;
    count = 1;
    while count<2%dtij > 1e-2      
        
        nnx = 2^count;
        
        nsubc = nnx^3;
        
        ddx = dx / nnx;
        ddy = dy / nnx;
        ddz = dz / nnx;
        
        zn = [Znode(K,I,J);Znode(K,I,J) - cumsum(ones(nnx,1)*ddz)];
        xn = [Xnode(K,I,J);Xnode(K,I,J) + cumsum(ones(nnx,1)*ddx)];
        yn = [Ynode(K,I,J);Ynode(K,I,J) + cumsum(ones(nnx,1)*ddy)];
        
        [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
        
        dobsx = xn - obsx ;
    
        dobsy = yn - obsy ;

        dobsz = obsz - zn ;
        
        [subcell,~,~] = topocheck(Xn,Yn,Zn,topo);
        
        index = find(subcell==1);
        
        % compute kernel for active cells
        [tx,ty,tz] = MAG3C_T_row(dobsx,dobsy,dobsz);
    
        tsuper(1,superC)        = tij(1,superC) + sum(tx(index));
        tsuper(1,superC+mcell)  = tij(1,superC+mcell) + sum(tx(index+nsubc));
        tsuper(1,superC+2*mcell)= tij(1,superC+2*mcell) + sum(tx(index+2*nsubc));
        tsuper(2,superC)        = tij(2,superC) + sum(ty(index));
        tsuper(2,superC+mcell)  = tij(2,superC+mcell) + sum(ty(index+nsubc));
        tsuper(2,superC+2*mcell)= tij(2,superC+2*mcell) + sum(ty(index+2*nsubc));
        tsuper(3,superC)        = tij(3,superC) + sum(tz(index));
        tsuper(3,superC+mcell)  = tij(3,superC+mcell) + sum(tz(index+nsubc));
        tsuper(3,superC+2*mcell)= tij(3,superC+2*mcell) + sum(tz(index+2*nsubc));
        
        dtij = abs(sum(sum(tsuper(:,nullcell==1)))-sum(sum(ttemp(:,nullcell==1))))/sum(sum(tsuper(:,nullcell==1)));
        
        ttemp = tsuper;
        
        count = count +1 ;

    end
    
    

end

Tx = tsuper(1,:);
Ty = tsuper(2,:);
Tz = tsuper(3,:);

