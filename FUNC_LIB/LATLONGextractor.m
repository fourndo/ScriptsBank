%Extract file name, LAT LON in each file
close all
clear all
%List all files in the folder and create structure

cd 'C:\PVK_Projects\3582_Rio_Tinto\data\Bingham 2008 2010';
listing=dir('C:\PVK_Projects\3582_Rio_Tinto\data\Bingham 2008 2010');


fileID = fopen('stations.txt','w');
for ii=3:size(listing,1);
    
%Take the name out of the structure listing(ii).name

newData1 = importdata(listing(ii).name);

% Create new variables in the base workspace from those fields.
vars = fieldnames(newData1);
for i = 1:length(vars)
    assignin('base', vars{i}, newData1.(vars{i}));
end

%Look in textdata for strings
fileid=listing(ii).name;
lon=textdata(strmatch('LON', textdata));
lat=textdata(strmatch('LAT', textdata));

fprintf(fileID,'%s %s %s\n',listing(ii).name,lon{1},lat{1});

end
fclose(fileID);



