function [G]=compG(X0,Z0,nX,nZ,dX,dZ,Rx_X,Rx_Z,Tx_X,Tx_Z)
%Create sensitivity matrix

%Pre-allocate space
G=sparse(1,nX * nZ + 1);
G(1)=1;

%Compute total ray path length
rangeX = Rx_X - Tx_X + 1e-11;
rangeZ = Rx_Z - Tx_Z + 1e-11;

%Get sign of increment (+ or -)
dirX = round( rangeX / abs(rangeX) );
dirZ = round( rangeZ / abs(rangeZ) );

%Compute incident angle
dip = dirZ * atan(abs(rangeZ/rangeX));

Lrange=(rangeX ^2 + rangeZ ^2) ^ 0.5;

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

if dirX > 0 %Ray travels East
    
   while round(dL*1e+5) < round(Lrange*1e+5) 
        if dirZ > 0 %Ray travels east-up
            
            if cellz == Tx_Z
                countZ=countZ - 1;
                cellz = cellz + dZ(countZ);
            end
            
            dx = cellx + dX(countX) - xin;
            dz = cellz - zin;
            
            %Find shortest path out of cell
            dlx= dx / cos(dip) ;
            dlz= dz / sin(dip);
            
            if dlx < dlz
            dl = dlx;
            G((countX-1)*nZ+countZ + 1) = dl;
            dL= dL + dl;
            
 
            else
                
            dl = dlz;
            G((countX-1)*nZ + countZ + 1) = dl;
            dL= dL + dl;
            end
            


        xin= xin + dl * cos(dip);
        zin= zin + dl * sin(dip);
        
            %Increment for next cell
            if xin >= cellx + dX(countX)
            cellx = cellx + dX(countX);
            countX = countX+1;
            end


            if dirZ > 0 && zin >= cellz
                countZ = countZ - 1;
                cellz = cellz + dZ(countZ);
                    
            elseif dirZ < 0 && zin <= cellz - dZ(countZ)
                countZ = countZ + 1;
                cellz = cellz - dZ(countZ);
                
            end 
            

        else %Ray travels east-down

            dx = cellx + dX(countX) - xin;
            dz = cellz - dZ(countZ) - zin;
            
            %Find shortest path out of cell
            dlx= dx / cos(dip) ;
            dlz= dz / sin(dip);
            
            if dlx < dlz
            dl = dlx;
            G((countX-1)*nZ+countZ + 1) = dl;
            dL= dL + dl;
            
 
            else
                
            dl = dlz;
            G((countX-1)*nZ + countZ + 1) = dl;
            dL= dL + dl;
            end
            


        xin= xin + dl * cos(dip);
        zin= zin + dl * sin(dip);
        
            %Increment for next cell
            if xin >= cellx + dX(countX)
            cellx = cellx + dX(countX);
            countX = countX+1;
            end


            if dirZ > 0 && zin >= cellz
                cellz = cellz + dZ(countZ);
                countZ = countZ - 1;    
            elseif dirZ < 0 && zin <= cellz - dZ(countZ)
                cellz = cellz - dZ(countZ);
                countZ = countZ + 1;
            end 
            

        end



            
   end
   

end

    
