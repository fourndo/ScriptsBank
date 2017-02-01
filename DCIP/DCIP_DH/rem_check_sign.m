% function rem_check_sign(check_file,obs_file)
%
% Load a UBC-DCIP3D "check_sign.txt" file and remove bad lines from the
% corresponding obs file

% Written by: D. Fournier
% Last update: January 3th, 2014

clear all
close all

%% USER INPUT
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Desktop\Stan';
check_file = [work_dir '\check_sign.txt'];
obs_file = [work_dir '\Arrow_PD_edit_dc.txt'];

% Name for edited output observation file
obs_file_out = [work_dir '\Arrow_PD_edit_dc_SignChecked.txt'];


flag = '639';
dtype = 4; % Number of columns (surface format [Ax, Bx, Mx, Nx] =4, general [Ax, Az, ..., Nx, Nz] = 8) 


%% ||| SCRIPT STARTS HERE |||
fsign = fopen(check_file,'r');
% Import check sign data
line = fgets(fsign);


% First fix the sign file, usually the first and second coordinates are
% stuck to each other
% Pre-allocate 2D array
checksign = zeros(1,dtype+2);
count = 0;

while line~=-1
    
    line = fgets(fsign);
    line = regexprep(line,['(' flag ')'],' $1');
    if isempty(str2num(line))==1
   
        continue     
    else
        
        count= count+1;
        
        
  
        line = fgets(fsign);
        line = regexprep(line,['(\.+[0-9]*)(' flag ')([0-9]*+\.)'],'$1 $2$3');
        
        data = str2num(line);
        
        % Extract data and round to precision
        data(1:dtype) = fix(data(1:dtype));
        checksign(count,1:dtype+1) = data;

        
    end
    
    
    line = fgets(fsign);
    
end

fclose(fsign);

% Open obsfile and remove lines with check sign flag
obs_in = fopen(obs_file,'r');
line = fgets(obs_in);

% Write new obs file
obs_out = fopen(obs_file_out,'w');
count =1;
while line~=-1
    
    data = str2num(line);
        
    if isempty(data)==0 && length(data)==dtype+2
        
        data(1:dtype) = fix(data(1:dtype));
        
%         X = fix(data(1));
%         Y = fix(data(2));
%         Z = fix(data(3));
%         
%         datum = fix(data(end-1) * 1e+4) / 1e+4;
        
        index = ones(size(checksign,1),1);
        for ii = 1:dtype
            index = (index.*(data(ii) == checksign(:,ii)))==1;
            if sum(index)==0
                break
            end
        end

        % If current line has a checksign flag, skip to next line
        if (sum(index)~=0)% && sum(Z == checksign(index,3))~=0 && sum(X == checksign(index,1))~=0 && sum(Y == checksign(index,2))~=0)
            
            % flip the sign
            checksign(index,end) = checksign(index,end)+1;
            data(end-1) = data(end-1)*-1;
            
            fprintf(obs_out,'%.2f ',data(1:end-2));

            fprintf(obs_out,'%12.8e %12.8e\n',data(end-1:end));

            line = fgets(obs_in);
            
            continue
            
        % Else continue
        else
            
           fprintf(obs_out,'%s',line);
           
        end
        
        
    else
        
        fprintf(obs_out,'%s',line);
        
    end
    
    line = fgets(obs_in);
    
    count = count+1;
    
    
end

fclose(obs_in);
fclose(obs_out);

