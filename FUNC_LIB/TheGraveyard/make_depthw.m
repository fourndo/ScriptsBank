function make_depthw(outdir,meshfile,mfile,layers)
% Create depth weighting model for DCIP
% Inputs:
% outdir: Directory to store the w.dat file
% meshfile: mesh used for the inversion
% mfile: model from the inversion. Used to get the discretized topography
% layers: number of depth weighting layers ( i.e. [64 32 16 8 4 2] )
%           size of the mesh must be at least equal to the number of layers
%
% Output:
% w.dat: weighting files in UBC format (Z,X,Y)

fprintf('Creating depth weighting matrix w.dat\n')
fprintf('Number of input layers: %i\n',length(layers));
ndv = 1e-8;

mesh = get_UBC_mesh(meshfile);

nx = mesh(1,1);
ny = mesh(1,2);
nz = mesh(1,3);
mcell = nx * ny * nz;

m = load(mfile);
% Pre-allocate space for weighting matrix

ws = ones(mcell,1);
wx = ones((nx-1)*ny*nz,1);
wy = ones(nx*(ny-1)*nz,1);
wz = ones(nx*ny*(nz-1),1);

countx = 1;
county = 1;
countm = 1;
counter = 1;
for jj = 1 : ny
    
    for ii = 1 : nx
                
        for kk = 1 : nz
                       
                
            if m(countm) ~= ndv
                
                if counter <= length(layers) 
                    
                    if ii < nx

                       wx(countx) = layers(counter);
                                              
                    end
                    
                    if jj < ny

                        wy(county) = layers(counter);

                    end
                    
                    counter = counter+1;

                end
               

                
            end
                
            if ii < nx
                
              countx = countx + 1;  
              
            end
            
            if jj < ny
                
              county = county + 1;  
              
            end 
            
            countm = countm + 1;
            
        end
        
        % Re-initilize for next column
        counter = 1;
        
    end
    
end
W = [ws;wx;wy;wz];
save([outdir '\w.dat'],'-ascii','W');

fprintf('Depth weigthting matrix created\n')
fprintf('--> %s\n',outdir)
