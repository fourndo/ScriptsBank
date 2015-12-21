function [G,xin,yin,zin]=getG_3D(Rx,phi,theta,x0,y0,z0,dX,dY,dZ,ZZ)

% Inputs
% Rx: [x,y,z] coordinates of receiver station
% phi: dipp angle
% theta: azimuth
% dX, dY, dZ: Mesh cell size (vectors)
% ZZ: Discretized topography (2D array)

nX = length(dX);
nY = length(dY);
nZ = length(dZ);

G = zeros( 1 , nX * nY * nZ );

%Initiate variables for cell steps
cellx=x0;
celly=y0;
cellz=z0;

%Search for starting cell position in grid (Top-West corner)
    countX=1;
    while cellx <= Rx(1)
        cellx = cellx + dX(countX);
        countX=countX + 1;
    end
    %Come back one cell
    countX=countX - 1;
    cellx = cellx - dX(countX);
    
    countY=1;
    while celly <= Rx(2)      
        celly = celly + dY(countY);
        countY= countY + 1;
    end
    countY=countY-1;
    celly = celly - dY(countY);
    
    countZ=1;
    while cellz >= Rx(3)      
        cellz = cellz - dZ(countZ);
        countZ= countZ + 1;
    end
    countZ=countZ-1;
    cellz = cellz + dZ(countZ);
    
    
%Compute first entry for remainder cell start
xin= Rx(1);
yin= Rx(2);
zin= Rx(3);

% Get direction of increment
if theta<0
    Xdir = -1;
else
    Xdir = 1;
end

if cos(theta)<0
    Ydir = -1;
else
    Ydir = 1;
end

% Compute until ray passes topography   
while zin < ZZ(countX,countY) 
%         if dirX > 0 %Ray travels east
            
            if cellx == xin && Xdir < 0
                countX=countX - 1;
                cellx = cellx - dX(countX);
            end
            
            if celly == yin && Ydir < 0
                countY=countY - 1;
                celly = celly - dY(countY);
            end
            
            if Xdir >= 0
                
            dx = (cellx + dX(countX)) - xin ;
            
            else
                
            dx = xin - cellx ;
            
            end
            
            if Ydir >= 0
                
            dy = (celly + dY(countY)) - yin ;
            
            else
                
            dy = yin - celly ;
            
            end
            
            dz = cellz - zin;
            
            %Find shortest path horizontally out
            dlx= abs( dx / sin(theta) ) ;
            
            dly= abs( dy / cos(theta) ) ;
            
            dlh = min([dlx dly]);
            
            %Find the shortest path out vertically
            dlz= abs( dz / cos(phi) ) ;
            
            dlhz = abs( dlh / sin(phi) ) ;
            
            % Compute overall shortest path 
            dl = min ([dlz dlhz]);
            
            % Update forward operator
            
            G((nZ) * (nX) * (countY-1) + nZ * (countX-1) + countZ) = dl;            

            xin= xin + Xdir * abs(dl * sin(theta) * sin(phi));
            yin= yin + Ydir * abs(dl * cos(theta) * sin(phi));
            zin= zin + dl * cos(phi);
        
            if xin > (x0 + sum(dX)) || xin < x0 || yin > (y0+sum(dY)) || yin < y0
                return
            end
            
            %Increment for next cell
            if xin >= cellx + dX(countX) && Xdir>0
            cellx = cellx + dX(countX);
            countX = countX+1;
            end

            if yin >= celly + dY(countY) && Ydir>0
            celly = celly + dY(countY);
            countY = countY+1;
            end

            if zin >= cellz && countZ>1
                
                countZ = countZ - 1;
                cellz = cellz + dZ(countZ);
                
            end 
            

end



            
