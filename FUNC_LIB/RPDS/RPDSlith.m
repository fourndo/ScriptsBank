%Extract file name, LAT LON in each file
close all
clear all
%List all files in the folder and create structure


load RPDS_Abbreviations
load RPDS_Era
load RPDS_Epoch
load RPDS_Period
load RPDS_Age
% load RPDS_Lithology
load RPDS_Master_Litho
% load RPDS_Common_rock
load RPDS_Formation
load RPDS_Atlas

% Import the file
file_in = importdata('litho_to_sort.csv');
rawData= file_in.textdata;

% Add character at the end of string for sorting
for ii=1:length(rawData)
    rawData{ii} = ['* ' rawData{ii} ' *'];
end
 

%% Find and replace common abbreviations
for nn=1:size(rawData,1)
    proData{nn,1}=rawData{nn,1};
    for kk=1:size(Abbreviation,1)
        proData{nn,1} = strrep(proData{nn,1}, Abbreviation{kk,1}, Abbreviation{kk,2});
    end
end

%% Extract all words from RPDS_Atlas
for ii=1:length(RPDS_Atlas)
match(:,ii)=regexpi(proData,RPDS_Atlas{ii});
end
data_filtered{size(match,1),1}='';

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       data_filtered{jj}=[data_filtered{jj} ' ' RPDS_Atlas{oo}];

       end
   end
   
end

clear match

% Change \s for space ' '
fileID = fopen('Data_filtered.txt','w');
for ii=1:length(data_filtered)
    if isempty(data_filtered{ii})==0
    data_filtered{ii}=strrep(data_filtered{ii}, '\s', '');
    end
   fprintf(fileID,'%s',data_filtered{ii});
   fprintf(fileID,'\n');
end
fclose(fileID);
%% Find and create master lithology list

for ii=1:length(Master_litho)
match(:,ii)=regexpi(data_filtered,Master_litho{ii});
end


fileID = fopen('Master_litho_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Master_litho{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

% %% Find and create lithology list
% for ii=1:length(Lithology)
% match(:,ii)=regexpi(data_filtered,Master_litho{ii});
% end
% 
% 
% fileID = fopen('Lithology_out.txt','w');
% 
% for jj=1:size(match,1);
% 
%    for oo=1:size(match,2);
%        
%        if isempty(match{jj,oo})==0
%        fprintf(fileID,'%s ',Lithology{oo});
%        end
%        
%    end
%    fprintf(fileID,'\n');
% end
% fclose(fileID);
% clear match
% 
% fid=fopen('RPDS_Atlas.dat','w');
% for jj=1:length(atlas)
%     fprintf(fid,'%s\n',atlas{jj});
% end
% fclose(fid)

% %% Find and create common rock list
% for ii=1:length(Common_rock)
% match(:,ii)=regexpi(proData,Common_rock{ii});
% end
% 
% 
% fileID = fopen('Rocks_out.txt','w');
% 
% for jj=1:size(match,1);
%     
%    for oo=1:size(match,2);
%        
%        if isempty(match{jj,oo})==0
%        fprintf(fileID,'%s ',Common_rock{oo});
%        end
%    end
%    fprintf(fileID,'\n');
% end
% fclose(fileID);
% clear match

%% Find and create Era
for ii=1:length(Era)
match(:,ii)=regexpi(data_filtered,Era{ii});
end


fileID = fopen('Era_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Era{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

%% Find and create Period
for ii=1:length(Period)
match(:,ii)=regexpi(data_filtered,Period{ii});
end


fileID = fopen('Period_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Period{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

%% Find and create Epoch
for ii=1:length(Epoch)
match(:,ii)=regexpi(data_filtered,Epoch{ii});
end


fileID = fopen('Epoch_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Epoch{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

%% Find and create Age
for ii=1:length(Age)
match(:,ii)=regexpi(data_filtered,Age{ii});
end


fileID = fopen('Age_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Age{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

%% Find and create Formation
for ii=1:length(Formation)
match(:,ii)=regexpi(data_filtered,Formation{ii});
end


fileID = fopen('Formation_out.txt','w');

for jj=1:size(match,1);
    
   for oo=1:size(match,2);
       
       if isempty(match{jj,oo})==0
       fprintf(fileID,'%s ',Formation{oo});
       end
   end
   fprintf(fileID,'\n');
end
fclose(fileID);
clear match

%% Create NO vector

no_vec=regexpi(proData,' no ');
no_vec=[no_vec regexpi(proData,' not ')];
no_vec=[no_vec regexpi(proData,' nor ')];

fileID = fopen('No_flag.txt','w');
for ii=1:size(no_vec,1)
   fprintf(fileID,'%d %d',no_vec{ii,1},no_vec{ii,2},no_vec{ii,3});
   fprintf(fileID,'\n');
end
fclose(fileID);