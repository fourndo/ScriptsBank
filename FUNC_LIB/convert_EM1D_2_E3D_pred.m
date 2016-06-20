function dmat = convert_EM1D_2_E3D_pred(obsfile,tc)
% function convert_dbs_2_mat(data,axis)
% Convert a E3D-TD file format to a 2D matrix. Only reformat the data for a
% single component at the time
%
% INPUT:
% obsfile   : EM1DTM obsfile file
% tc        : Time channels expected by E3D file
% outfile   : Name of the outputed predicted E3D file
%
% OUTPUT:
% outfile

% Load EM1D obs file
data = read_EM1DTM_obs(obsfile);

dmat = NaN(size(data{5},1)*size(data{5}{1},1)*length(tc),13);


count = 0;
% Loop through all the stations and write file
for ii = 1 : size(data{5},1)
    
    for jj = 1 : size(data{5}{ii},1)
        
       tin = data{5}{ii}{jj}{5}(:,1);
       
       % Write to file
       for kk = 1 : length(tc)
           count = count + 1;
           index = tc(kk) == tin;
           
               
               dmat(count,1:4) = [data{5}{ii}{jj}{2}(1),data{5}{ii}{jj}{2}(2),...
                   data{5}{ii}{jj}{2}(3), tc(kk)];
               
           if sum(index) == 1
               dmat(count,end) = data{5}{ii}{jj}{7}(index,1); 
           end

            
       end
       
       
    end
    
end

    
    
