function [H, BI, BD, MI, MD, obsx, obsy, obsz, d, wd] = read_MAG3D_obs(Obsfile)
% [ObsX, ObsY, ObsZ, d, Wd, H, I, D] = read_UBCobs(Obsfile)
% Load UBC observations and generate error weighting matrix

fid = fopen(Obsfile,'r');

% Extract Incl, Declin and field strength
line = fgets(fid); 
data = regexp(strtrim(line),'\s*','split');
BI = str2double(data{1});
BD = str2double(data{2});
H = str2double(data{3});

% Adjust D:azimuth from North to Cartesian
% [D] = mod(450-BD,360);


% Extract Incl, Declin of magnetization and type of data
line = fgets(fid);
data = regexp(strtrim(line),'\s*','split');
MI = str2double(data{1});
MD = str2double(data{2});
dtype = str2double(data{3}); 

% Extract number of data for each components
line = fgets(fid);
data = regexp(strtrim(line),'\s*','split');


if dtype==1 
    ndata = str2double(data{1});
    
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
    
elseif dtype==0 && mod(1,str2double(data{1})/3)==1

    ndata = str2double(data{1})/3;
    fprintf('It appears that the data has 3 components\n');
    fprintf('Program will assume ndata/3 \n');
    Bx = zeros( ndata , 1 ) ;
    wdx = zeros( ndata , 1 ) ;

    By = zeros( ndata , 1 ) ;
    wdy = zeros( ndata , 1 ) ;

    Bz = zeros( ndata , 1 ) ;
    wdz = zeros( ndata , 1 ) ;

    obsx = zeros( ndata , 1 ) ;
    obsy = zeros( ndata , 1 ) ;
    obsz = zeros( ndata , 1 ) ;

    line = fgets(fid);
    data = str2num(line);

    for ii = 1: ndata

        obsx(ii) = data(1);
        obsy(ii) = data(2);
        obsz(ii) = data(3);

         if length(data) > 5

            Bx(ii) = data(6);

            if length(data) > 6

                wdx(ii) = (abs(data(7)));

            else

                wdx(ii) = 1;

            end

        else

            Bx(ii) = 0;

        end

        line = fgets(fid);
        data = str2num(line);

    end

    % Import y-component
    for ii = 1: ndata

        if isempty(data)==1

            fprintf('It appears that the Y component is not in the file. Please revise!\n');

        end

         if length(data) > 5

            By(ii) = data(6);

            if length(data) > 6

                wdy(ii) = (abs(data(7)));

            else

                wdy(ii) = 1;

            end

        else

            By(ii) = 0;

        end

        line = fgets(fid);
        data = str2num(line);

    end

    % Import z-component
    for ii = 1: ndata

        if isempty(data)==1

            fprintf('It appears that the Z component is not in the file. Please revise!\n');

        end

        if length(data) > 5

            Bz(ii) = data(6);

            if length(data) > 6

                wdz(ii) = (abs(data(7)));

            else

                wdz(ii) = 1;

            end

        else

            Bz(ii) = 0;

        end

        if ii~=ndata
        line = fgets(fid);
        data = str2num(line);
        end

    end

    d   = [Bx;By;Bz];
    wd  = [wdx;wdy;wdz];

    fclose(fid);
end   

fclose all;