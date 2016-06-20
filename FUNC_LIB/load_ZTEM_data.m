function [data,ndv,dtype] = load_ZTEM_data(obsfile)
% Function loads a UBC-ZTEM file and output a data matrix
%
% INPUT:
% Obsfile: Location and name of file
%
% OUTPUT
% data: [Freq, x, y, z, Xreal, std, Ximag, std, Yreal, std, Yimag, std]

data = [];
% Open file
fid = fopen(obsfile,'r');

% Get data type
line = fgets(fid);
arg = regexp(line,'\s','split');
dtype = arg{2};

% No-data-value
line = fgets(fid);
arg = regexp(line,'\s','split');
ndv = arg{2};

% Skip to beginning of data
line = fgets(fid);

while line ~= -1
    
    if isempty(strtrim(line))
        line = fgets(fid);
        continue
    end
    
    % Get frequency
    arg = regexp(line,'\s','split');
    freq = str2num(arg{1});
    
    % Get number of lines
    line = fgets(fid);
    nrow = str2num(line);
    
    line = fgets(fid);
    d = zeros(nrow,12);
    
    for ii = 1 : nrow
        
        d(ii,:) = [freq str2num(line)];
        
        line = fgets(fid);
        
    end
    
    data = [data;d];
    
end
