function [G]=getG_2D(Rx,phi,Xo,Zo,dX,dZ,ZZ)

% Inputs
% Rx: [x,y,z] coordinates of receiver station
% phi: dipp angle
% theta: azimuth
% dX, dY, dZ: Mesh cell size (vectors)
% ZZ: Discretized topography (2D array)

nX = length(dX);

nZ = length(dZ);

G = zeros( 1 , nX * nZ );

%Initiate variables for cell steps
cellx=Xo;

cellz=Zo;

%Search for starting cell position in grid (Top-West corner)
    countX=1;
    while cellx <= Rx(1)
        cellx = cellx + dX(countX);
        countX=countX + 1;
    end
    %Come back one cell
    countX=countX - 1;
    cellx = cellx - dX(countX);
    
    
    countZ=1;
    while cellz > Rx(2)      
        cellz = cellz - dZ(countZ);
        countZ= countZ + 1;
    end
    countZ=countZ-1;
    cellz = cellz + dZ(countZ);
    
    
%Compute first entry for remainder cell start
xin= Rx(1);
zin= Rx(2);

% Get direction of increment
if (phi)<0
    Xdir = -1;
else
    Xdir = 1;
end

% Compute until ray passes topography   
while zin < ZZ(countX) 
%         if dirX > 0 %Ray travels east
            
            if cellx == xin && Xdir < 0
                countX=countX - 1;
                cellx = cellx - dX(countX);
            end
            
            
            if Xdir >= 0
                
            dx = (cellx + dX(countX)) - xin ;
            
            else
                
            dx = xin - cellx ;
            
            end
            
            
            dz = cellz - zin;
            
            %Find shortest path out of cell
            dlx= abs( dx / sin(phi) ) ;
            dlz= abs( dz / cos(phi) ) ;
            
            % Compute overall shortest path 
            dl = min ([dlx dlz]);
            
            % Update forward operator
            G((countX-1)*nZ + countZ) = dl;
            
            xin= xin + dl * sin(phi);
            zin= zin + dl * cos(phi);

        
            if xin > (Xo + sum(dX)) || xin < Xo
                return
            end
            
            %Increment for next cell
            if xin >= cellx + dX(countX) && Xdir>0
            cellx = cellx + dX(countX);
            countX = countX+1;
            end

            if zin >= cellz && countZ>1
              
                
                countZ = countZ - 1;
                cellz = cellz + dZ(countZ);
                
            end 
            

end



            
