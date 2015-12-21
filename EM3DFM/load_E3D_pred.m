function [d_pred] = load_E3D_pred(predfile)

fid = fopen(predfile,'r');

% ndata = size(dwndata.data,1);
line = fgets(fid);
d_pred = zeros(1,16);
freq = [];
count = 1;
while line~=-1
    
    if isempty(regexp(line,'frequency','match'))==0
        
        temp = regexp(line,'[=]','split');
        freq = str2double(temp{2});
        
    end
        
    if isempty(regexp(line,'%','match'))==1 && isempty(regexp(line,'\d\s','match'))==0
        
        d_pred(count,1) = freq;
        d_pred(count,2:end) = str2num(line);
        count = count+1;
        
    end
    
    line = fgets(fid);
    
    
end

fclose(fid);
