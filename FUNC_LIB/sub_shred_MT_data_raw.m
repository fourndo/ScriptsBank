function[num_freq_all_sites, data, stn_name_list]=sub_shred_MT_data(file)
%% Function reformats MT data from blocks to full column ascii

% site1
% freq rot Zxxr Zxxi  ZxxVar  Zxyr .. .. Zyxr .. .. Zyyr .. ..
% : :
% : :
% : :
% site2
% : :


%%
format long e

max_num_lines=30000;   %some large number that we are sure is greater than the number of lines in the file

data_block=6;   % number of column in each data block

%% counts lines and number of sites per line

% fid=fopen(file,'rt');
% match_freq=('Frequency');
% d=fscanf(fid,'%s');
% a=findstr(d, match_freq);
% num_sites=length(a); %#ok<NOPTS>
% 
% fclose(fid);

%% reads and parses data
nb_files=20;
file_list=ls('C:\Projects\3893_Barrick_Goldstrike\Location_Data\Goldstrike\Shredit');
cd Shredit
fid2=fopen('3893_Goldstrike.dat','w');

nb_files=size(file_list,1)-2;

for oo=1:nb_files
   AA= strtok(file_list(oo+2,:),' ');
fid=fopen(AA,'rt');

line_num=AA(4:5);

num_freq_all_sites=[];
data=[];
stn_name_list={};

match_freq=('Frequency');
match_freq2=('Freq~i');
match_site=('SITE');
count=0;

count_stations=0;
while 1
    line=fgets(fid); %gets next line
    if line==-1
        break
    else

        count=count+1;
        fprintf('read line %i\n',count);

        if length(line) >= 5
            test_freq_match=strcmp(line(1:length(match_freq)),match_freq);
            test_freq_match2=strcmp(line(1:length(match_freq2)),match_freq2);
            test_site_match=strcmp(line(1:length(match_site)),match_site);
        else
            test_freq_match=0;
            test_freq_match2=0;
            test_site_match=0;
        end


        if test_site_match==1

            [str,rem]=strtok(line,':');
            rem=deblank(rem);
            stn_name=rem(2:end);
           % keyboard;
            
        end

        if test_freq_match==1 | test_freq_match2==1
            
            count_stations=count_stations+1;
            % get number of frequencies
            num_freq_str=fgetl(fid);
            count=count+1;
            fprintf('read line %i\n:',count);
            [name num_freq]=strread(num_freq_str, '%s %d,');
            num_freq_all_sites=[num_freq_all_sites; num_freq];
            num_freq_lines=ceil(num_freq/data_block);   %number of data lines to read for each block

            %get site data
            block=[];
            for k=1:14  %number of blocks per site for this type of data
                subblock=[];  j=[];
                for j=1:num_freq_lines
                    subblock_temp=str2num(fgetl(fid))';
                    count=count+1;
                    fprintf('read line %i\n:',count);
                    subblock=[subblock; subblock_temp];
                end

                block=[block subblock];
                temp=fgets(fid);
                count=count+1;
                fprintf('read line %i\n:',count);
                temp=fgets(fid); % move to top of next block
                count=count+1;
                fprintf('read line %i\n:',count);
            end


            %keyboard;

            all_data=block;%

            %Convert units feild to SI
            
            all_data(:,3:end)=all_data(:,3:end)/(10000/(4*pi));
            
            [A,stn_num]=strread(stn_name,'%s%d','delimiter','_');
           
            

            
            for jj=1:size(all_data,1)
                fprintf(fid2,'%s %d %12.10e %12.10e %12.10e %12.10e %12.10e %12.10e %12.10e %12.10e %12.10e\n',...
                    line_num,stn_num,all_data(jj,1),all_data(jj,3),all_data(jj,4),all_data(jj,6),all_data(jj,7),...
                    all_data(jj,9),all_data(jj,10),all_data(jj,12),all_data(jj,13))
            end
            
            % Apparent resistivity and phase
            mu_0=4*pi*1e-7;

%             all_data(:,15)=( (all_data(:,3).^2) + (all_data(:,4).^2) ) ./ (all_data(:,1)*mu_0*2*pi);   %appres [Ohm-m] xx
%             all_data(:,17)=( (all_data(:,6).^2) + (all_data(:,7).^2) ) ./ (all_data(:,1)*mu_0*2*pi);   %appres [Ohm-m] xy
%             all_data(:,19)=( (all_data(:,9).^2) + (all_data(:,10).^2) ) ./ (all_data(:,1)*mu_0*2*pi);   %appres [Ohm-m] yx
%             all_data(:,21)=( (all_data(:,12).^2) + (all_data(:,13).^2) ) ./ (all_data(:,1)*mu_0*2*pi);   %appres [Ohm-m] yy
% 
%             all_data(:,16)= (180/pi)*(atan( all_data(:,4)./all_data(:,3))) ;   %phase [deg] xx
%             all_data(:,18)= (180/pi)*(atan( all_data(:,7)./all_data(:,6))) ;   %phase [deg] xy
%             all_data(:,20)= (180/pi)*(atan( all_data(:,10)./all_data(:,9))) ;   %phase [deg] yx
%             all_data(:,22)= (180/pi)*(atan( all_data(:,13)./all_data(:,12))) ;   %phase [deg] yy

%             plotcurves(all_data,stn_name)

            %pause

            S=size(block);
            for ii=1:1:S(1)
                stn_name_list=[stn_name_list, stn_name];
            end

            data=[data; block];

        end
    end
end
end
fclose(fid);
fclose(fid2);


%1 Frequency
%2 Rotation
%3 Zxxr
%4 Zxxi
%5 Zxx variance
%6 Zxyr
%7 Zxyi
%8 Zxy variance
%9 Zyxr
%10 Zyxi
%11 Zyx variance
%12 Zyyr
%13 Zyyi
%14 Zyy variance
%15 Rhoxx
%16 Pxx
%17 Rhoxy
%18 Pxy
%19 Rhoyx
%20 Pyx
%21 Rhoyy
%22 Pyy





% checks
% num_freq_all_sites
% sum(num_freq_all_sites)
% size(data)

%% save
% save freq.dat num_freq_all_sites -ascii
% save mt.dat data -ascii 