function [mesh]=get_UBC_mesh(meshfile)
% Read UBC mesh file and extract parameters
% Works for the condenced version (20 * 3) --> [20 20 20] 
fid=fopen(meshfile,'rt');


% Go through the log file and extract data and the last achieved misfit
for ii=1:5         	
line=fgets(fid);
    
    % First line: number of cells in i, j, k 
    if ii==1
        mesh(1,:) = str2num(line);
    
    
    % Second line: origin coordinate (X,Y,Z)
    elseif ii==2
        mesh(2,:) = str2num(line);
    
    % Other lines for the dX, dY ,dZ
    else
        
        % Split the line
        var = regexp(strtrim(line),'\s*','split');
        vec = zeros(1,mesh(1,ii-2));
        
        count = 1;
        for jj = 1 : length(var)

            if isempty(regexp(var{jj},'[*]','match'))==1

                vec(count) = str2double(var{jj});
                count = count + 1;
            else

                temp = regexp(var{jj},'*','split');
                dl = str2double(temp(2));
                ndx = str2double(temp(1));

                vec(count:count+ndx-1) = ones(1,ndx) * dl;
                count = count+ndx;
            end

        end
        
        
            
        mesh(ii,1:mesh(1,ii-2)) = vec;
            
         
    end
    
      
end

fclose(fid);