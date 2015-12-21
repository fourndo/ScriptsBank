function [Wr] = get_Wr(mesh,data,topo_model)
%Create forward operator

dX = mesh(3,mesh(3,:)~=0)';
dY = mesh(4,mesh(4,:)~=0)';
dZ = mesh(5,mesh(5,:)~=0)';

nX = mesh(1,1);
nY = mesh(1,2); 
nZ = mesh(1,3);

X0 = mesh(2,1);
Y0 = mesh(2,2); 
Z0 = mesh(2,3);

R0 = min( [min(dX) min(dY) min(dZ)] ) / 4;
mcell = nX * nY * nZ;

%Pre-allocate for weighting matrix
wr=zeros(mcell,1);
Wr=zeros(mcell,1);
V=zeros(mcell,1);


for oo = 1:size(data,1)
    
    ObsX=data(oo,1);
    ObsY=data(oo,2);
    ObsZ=data(oo,3);
    
    count = 1;
    for jj = 1 : nY

        %First compute de location of the center of the cell
        Y = Y0 + sum(dY(1:jj)) - dY(jj) /2;

            for ii = 1 : nX

                X = X0 + sum(dX(1:ii)) - dX(ii) /2;

                    for kk = 1: nZ

                        Z = Z0 - sum(dZ(1:kk)) + dZ(kk) /2;

            %             G(count)=-G(count);
                        R= ((ObsX - X) ^ 2 + (ObsY - Y)^2 + (ObsZ - Z)^2)^0.5;

                        %Compute volume of cell
                        V(count) = dX(ii)*dY(jj)*dZ(kk);

                        % Compute the distance weighting
                        wr(count) = ( V(count) / (R + R0) ^ 3 ) ^ 2 ;

                        % Compute the forward operator

                        count = count + 1;

                   end

           end

    end
    
    Wr = Wr + wr;
    
    if mod(oo,100)==0
        oo
    end
    
end

Wr=Wr.^(1/4);

% Normalize depth weighting with the largest value

Wr = Wr./sqrt(V);

Wr = Wr .* topo_model;

Wr = Wr./(max(Wr));
