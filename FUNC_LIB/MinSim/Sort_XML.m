function [str,id] = Sort_XML(inputfile)
% Function reads in an XML file from GOcad and output a cell array of
% content.

fid = fopen(inputfile,'r');

line = fgets(fid);
count = 0;
while line ~=-1
    
    if isempty(regexp(line,'(<classification>)','match')) == 0
        
        line = fgets(fid);
        line = fgets(fid);
        line = fgets(fid);
        
    end
    
    if isempty(regexp(line,'(<name>)','match')) == 0
        
        temp = regexp(line,'>\w*<','match');
        str{count} = temp{1}(2:end-1);
        
    end    
        
    if isempty(regexp(line,'(<code>)','match')) == 0 
        count = count + 1;
        temp = regexp(line,'\d*','match');
        id(count) = str2num(temp{:});
        
    end
    
    line = fgets(fid);
    
end
        

