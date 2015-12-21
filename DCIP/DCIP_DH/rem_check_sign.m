% function rem_check_sign(check_file,obs_file)
%
% Load a UBC-DCIP3D "check_sign.txt" file and remove bad lines from the
% corresponding obs file

% Written by: D. Fournier
% Last update: January 3th, 2014

clear all
close all

%% USER INPUT
work_dir = 'C:\Projects\4160_Abitibi_Windfall\Inversion';
check_file = [work_dir '\check_sign.txt'];
obs_file = [work_dir '\Obs_AcrossHole_Tx1_DC_err.dat'];

% Name for edited output observation file
obs_file_out = [work_dir '\Obs_AcrossHole_Tx1_DC_err_EDT.dat'];


%% ||| SCRIPT STARTS HERE |||
fsign = fopen(check_file,'r');
% Import check sign data
line = fgets(fsign);

% Pre-allocate 2D array
checksign = zeros(1,8);
count = 0;
while line~=-1
    
    
    if isempty(str2num(line))==1
        
        
        line = fgets(fsign);
        line = regexprep(line,'\.[0-9]+(543)[0-9]+\.',' 543');
        
        continue
        
    else
        
        count= count+1;
        line = fgets(fsign);
        
        repstr = regexp(line,'543[0-9]+\.','match');
        repstr{1} =[' ' repstr{1}];
        repstr{2} =[' ' repstr{2}];
        
        line = regexprep(line,'543[0-9]+\.',repstr{1});
        
        data = str2num(line);
        
        data(5) = str2num(repstr{2});
        % Extract data and round to precision
%         data([1 3 4 6]) = data([1 3 4 6])*1e+1;
        data(end) = data(end)*1e+2;
        checksign(count,1:7) = fix(data);

        
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
        
    if isempty(data)==0 && length(data)==8
        
        data = data(1:end-1);
%         data([1 3 4 6]) = data([1 3 4 6])*1e+1;
        data(end) = data(end)*1e+2;
        data = fix(data);
        
%         X = fix(data(1));
%         Y = fix(data(2));
%         Z = fix(data(3));
%         
%         datum = fix(data(end-1) * 1e+4) / 1e+4;
        
        index = ((data(1) == checksign(:,1)) .* (data(2) == checksign(:,2)) .*...
             (data(3) == checksign(:,3)) .* (data(4) == checksign(:,4)) .*...
             (data(6) == checksign(:,6)) .* (data(5) == checksign(:,5)) .*...
             (data(7) == checksign(:,7)))==1;
        % If current line has a checksign flag, skip to next line
        if (sum(index)~=0)% && sum(Z == checksign(index,3))~=0 && sum(X == checksign(index,1))~=0 && sum(Y == checksign(index,2))~=0)
            
            
            checksign(index,8) = checksign(index,8)+1;
            
            
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

