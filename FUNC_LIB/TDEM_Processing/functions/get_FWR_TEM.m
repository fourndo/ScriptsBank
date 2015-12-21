function [FWR] = get_FWR_TEM(FWR_times,FWRfile,nstn)
% Extract dBdt from FWR file

FWR.XYZ = zeros(nstn,3);
FWR.dBx = zeros(nstn,length(FWR_times));
FWR.dBy = zeros(nstn,length(FWR_times));
FWR.dBz = zeros(nstn,length(FWR_times));


fid = fopen(FWRfile,'rt');
line=fgets(fid);
line=fgets(fid);
line=fgets(fid);
count = 1;
counter = 1;
while line~=-1
        
    temp = str2num(line);
    
    index = find(temp(4)==FWR_times);
    
    FWR.XYZ(count,1:3) = temp(1:3);
    
    FWR.dBx(count,index) = temp(11);
        
    FWR.dBy(count,index) = temp(12);
    
    FWR.dBz(count,index) = temp(13);
    
    counter = counter + 1;
    
    if counter > length(FWR_times)
        count = count +1;
        counter = 1;
    end
    
    line=fgets(fid);
    
end