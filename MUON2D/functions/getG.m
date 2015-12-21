function [G]=getG(X0,Z0,nX,nZ,dX,dZ,Tx_X,Tx_Z,dip)
% Inputs
% X0,Z0: Coordinate of top-west corner of the grid
% nX,nZ: Number of cells in each direction
% dX, dZ: Mesh cell size (vectors)
% Tx_X, Tx_Z: Coordinates of entry point for incomcing ray
% dip: angle of incidence (rad from vertical)

% function [G]=getG(Xmax,Zmax,Xo,Zo,dx,dz,dip)
%Create sensitivity matrix

% %Pre-allocate space
% G=zeros(1,(Xmax/dx)*(Zmax/dz));
% 
% %For now since all receivers are at same elevation
% Zo=Zmax;
% 
% %Create empty matrix for address
% model=zeros(Zmax/dz,Xmax/dx);
% 
% dl=dz/cosd(dip);
% 
% 
% while floor(Zo/dz)>0   
%     
%     Xf=Xo-dl*sind(dip);
%     Zf=Zo-dl*cosd(dip);
%     
%     if floor(Xf/dx)~=floor(Xo/dx)
%         if dip>0
%             dl_short=(Xo/dx-floor(Xo/dx))/sind(dip)*dx;
%             G(sub2ind(size(model),Zo/dz,ceil(Xo/dx)))=dl_short;
%             
%             dl_long=(floor(Xo/dx)-Xf/dx)/sind(dip)*dx;
%             G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl_long;
%         else
%             dl_short=abs((ceil(Xo/dx)-Xo/dx)/sind(dip))*dx;
%             G(sub2ind(size(model),Zo/dz,ceil(Xo/dx)))=dl_short;
%             
%             dl_long=abs((Xf/dx-ceil(Xo/dx))/sind(dip))*dx;
%             G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl_long;
%         end
%     else
%         G(sub2ind(size(model),Zo/dz,ceil(Xf/dx)))=dl;%*(0.5*cos(Zo*2*pi/Zmax+pi)+0.5);
%     end
%     Xo=Xf;
%     Zo=Zf;
% end
% clear model;
% end

G=zeros(1,nX * nZ);

%Compute total ray path length
% rangeX = Rx_X - Tx_X + 1e-11;
L = (Tx_Z - (Z0- sum(dZ))) / cos(dip);

% %Get sign of increment (+ or -)
Xdir = round( dip / abs(dip) );

%Initiate variables for cell steps
dL=0;
cellx=X0;
cellz=Z0;

%Search for starting cell position in grid (Top-West corner)
    countX=1;
    while cellx <= Tx_X
        cellx = cellx + dX(countX);
        countX=countX + 1;
    end
    %Come back one cell
    countX=countX - 1;
    cellx = cellx - dX(countX);
    
    
    countZ=1;
    while cellz >= Tx_Z      
        cellz = cellz - dZ(countZ);
        countZ= countZ + 1;
    end
    countZ=countZ-1;
    cellz = cellz + dZ(countZ);
    
    
%Compute first entry for remainder cell start
xin= Tx_X;
zin= Tx_Z;


% Compute until full ray path length = L (+- 1e-6 precision)   
while round(dL*1e+6) < round(L * 1e+6) 
%         if dirX > 0 %Ray travels east
            
            if cellx == xin && Xdir < 0
                countX=countX - 1;
                cellx = cellx - dX(countX);
            end
            
            if Xdir > 0
                
            dx = (cellx + dX(countX)) - xin ;
            
            else
                
            dx = xin - cellx ;
            
            end
            
            dz = zin - (cellz - dZ(countZ)) ;
            
            %Find shortest path out of cell
            dlx= abs( dx / sin(dip) ) ;
            dlz= abs( dz / cos(dip) ) ;
            
            if dlx < dlz
                
            dl = dlx;  
            
            else
                
            dl = dlz;
            
            end
            
            G((countX-1)*nZ + countZ) = dl;
            
            dL= dL + dl;

            xin= xin + dl * sin(dip);
            zin= zin - dl * cos(dip);
        
            if round(dL*1e+6) == round(L * 1e+6)
                return
            end
            
            %Increment for next cell
            if xin >= cellx + dX(countX) && Xdir>0
            cellx = cellx + dX(countX);
            countX = countX+1;
            end


            if zin <= cellz - dZ(countZ)
                
                countZ = countZ + 1;
                cellz = cellz - dZ(countZ);
                
            end 
            

%         else %Ray travels west
% 
%             dx = cellx + dX(countX) - xin;
%             dz = cellz - dZ(countZ) - zin;
%             
%             %Find shortest path out of cell
%             dlx= dx / cos(dip) ;
%             dlz= dz / sin(dip);
%             
%             if dlx < dlz
%             dl = dlx;
%             G((countX-1)*nZ+countZ + 1) = dl;
%             dL= dL + dl;
%             
%  
%             else
%                 
%             dl = dlz;
%             G((countX-1)*nZ + countZ + 1) = dl;
%             dL= dL + dl;
%             end
%             
% 
% 
%         xin= xin + dl * cos(dip);
%         zin= zin + dl * sin(dip);
%         
%             %Increment for next cell
%             if xin >= cellx + dX(countX)
%             cellx = cellx + dX(countX);
%             countX = countX+1;
%             end
% 
% 
%             if dirZ > 0 && zin >= cellz
%                 cellz = cellz + dZ(countZ);
%                 countZ = countZ - 1;    
%             elseif dirZ < 0 && zin <= cellz - dZ(countZ)
%                 cellz = cellz - dZ(countZ);
%                 countZ = countZ + 1;
%             end 
%             
% 
%         end
end



            
