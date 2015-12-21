function [G]=compG(X0,Z0,nX,nZ,dX,dZ,Rx_X,Rx_Z,Tx_X,Tx_Z)
%Create sensitivity matrix

%Pre-allocate space
G=zeros(1,nX * nZ);

%Compute total ray path length
rangeX = Rx_X(jj) - Tx_X(ii);
rangeZ = Rx_Z(jj) - Tx_Z(ii);

%Compute incident angle
angl_in = atan(rangeZ/rangeX);

Lrange=(rangeX ^2 + rangeZ ^2) ^ 0.5;

%Initiate variables for cell steps
dl=0;
count=1;
dL=0;
cellx=X0;
cellz=Z0;

%Search for starting cell position in grid
    countX=1;
    while cellx < Rx_X
        cellx = cellsx + dX(countX);
        countX=countX+1;
    end
    
    countZ=1;
    while cellz < Rx_Z
        cellz = cellz + dZ(countZ);
        countZ=countZ+1;
    end
        
while dL < Lrange   
    
    
    Xf=Xo-dl*sind(dip);
    Zf=Zo-dl*cosd(dip);
    
    if floor(Xf/dx)~=floor(Xo/dx)
        if dip>0
            dl_short=(Xo/dx-floor(Xo/dx))/sind(dip)*dx;
            G(sub2ind(size(model),Zo/dz,ceil(Xo/dx)))=dl_short;
            
            dl_long=(floor(Xo/dx)-Xf/dx)/sind(dip)*dx;
            G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl_long;
        else
            dl_short=abs((ceil(Xo/dx)-Xo/dx)/sind(dip))*dx;
            G(sub2ind(size(model),Zo/dz,ceil(Xo/dx)))=dl_short;
            
            dl_long=abs((Xf/dx-ceil(Xo/dx))/sind(dip))*dx;
            G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl_long;
        end
    else
        G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl;%*(0.5*cos(Zo*2*pi/Zmax+pi)+0.5);
    end
    Xo=Xf;
    Zo=Zf;
end
clear model;
end