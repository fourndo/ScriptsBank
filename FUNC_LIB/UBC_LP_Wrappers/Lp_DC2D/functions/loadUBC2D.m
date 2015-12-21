function [model] = loadUBC2D(model_file)
% Load a UBC-DCIP2D file and reformat in vector

fid=fopen(model_file,'rt');
line=str2num(fgets(fid));
nX = line(1);
nZ = line(2);
model = zeros(nZ,nX);
countx=0;
countz=1;
% Go through the log file and extract data and the last achieved misfit
while countz<=nZ
    line=fgets(fid); %gets next line 
    
    m_line = str2num(line);
    
    
    model(countz,:) = m_line;
    
    countz = countz+1;
    
end
        
model = reshape(model,nX*nZ,1);    
    
