function [H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_TMI(obsfile)
% [ObsX, ObsY, ObsZ, d, Wd, H, I, D] = read_UBCobs(Obsfile)
% Load UBC observations and generate error weighting matrix

fid = fopen(obsfile,'r');
line = fgets(fid);

line = regexp(line,'!!','split');

data = str2num(line{1});

I = data(1);
Dazm = data(2);
H = data(3);

% Adjust D:azimuth from North to Cartesian
[D] = mod(450-Dazm,360);

% Need to write code here to look for mode 2 flag and outpus components
line = fgets(fid);

line = fgets(fid);

line = regexp(line,'!!','split');
ndata = str2double(line{1});

d = zeros( ndata , 1 ) ;
wd = zeros( ndata , 1 ) ;
obsx = zeros( ndata , 1 ) ;
obsy = zeros( ndata , 1 ) ;
obsz = zeros( ndata , 1 ) ;

for ii = 1 : ndata
    
    line = fgets(fid);
    
    data = str2num(line);
    obsx(ii) = data(1);
    obsy(ii) = data(2);
    obsz(ii) = data(3);
    
    % If data is present
    if length(data)>3
        
        d(ii) = data(4);

        % If uncertainty is present
        if length(data)>4

            wd(ii) = (abs(data(5)));
        else
            
            wd(ii) = 1;
            
        end  
        
    else
        
        d(ii) = 0;
        
    end
    
end

fclose(fid);