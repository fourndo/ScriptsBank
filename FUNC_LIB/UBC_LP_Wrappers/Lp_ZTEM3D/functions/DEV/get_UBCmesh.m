function [mesh]=get_UBCmesh(meshfile)
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
        
        delm_entries = find(isspace(line)==1);
        
        dX = str2num(line(1:delm_entries(1)));
        
        for jj = 1:length(delm_entries)-1;
            
            read_str = line(delm_entries(jj):delm_entries(jj+1));
            
            index = strfind(read_str,'*');
            
            if isempty(index)==1
                
                dX =  [dX str2num(read_str)];
                
            else
                
                cell_size = str2num(read_str(index+1 : end));
                numcells = str2num(read_str(1 : index-1));
                
                dX =  [dX ones(1,numcells) * cell_size];
            end
            
        end
        
        mesh(ii,1:mesh(1,ii-2)) = dX;
         
    end
    
end