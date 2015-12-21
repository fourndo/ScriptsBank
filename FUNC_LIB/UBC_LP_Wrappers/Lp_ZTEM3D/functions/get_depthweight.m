function [wxdepth,wydepth] = get_depthweight(layers,nZ,nX,nY,topo_model)
% Create depth weighting function using a "layer" pattern.


% Apply depth weighting
wxdepth = ones(nZ,nX-1,nY);
wydepth = ones(nZ,nX,nY-1);
topo_model = reshape(topo_model,nZ,nX,nY);

    for jj = 1:nY
        for ii = 1:nX-1

                colm = topo_model(:,ii,jj) == 1;
                numcell = sum(colm);
                wxdepth(colm,ii,jj) = [layers;ones(numcell-length(layers),1)];                           

        end
    end

for jj = 1:nY-1
    for ii = 1:nX

            colm = topo_model(:,ii,jj) == 1;
            numcell = sum(colm);
            wydepth(colm,ii,jj) = [layers;ones(numcell-length(layers),1)];                           

    end
end

wxdepth = reshape(wxdepth,nZ*(nX-1)*nY,1);
wydepth = reshape(wxdepth,nZ*nX*(nY-1),1);