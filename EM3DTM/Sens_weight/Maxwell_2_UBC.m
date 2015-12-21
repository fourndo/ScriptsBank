%% Import xyz EM SQUID data exported from Maxwell and generate a UBC format datafile
%


%% Import DATA format
% %X  Y  LEVEL  ELEV  STATION  COMPONENT  CURRENT  DELAY  Ch1  Ch2  Ch3 Ch4  Ch5  Ch6  Ch7  Ch8  Ch9  Ch10  Ch11  Ch12  Ch13  Ch14  Ch15  Ch16  Ch17  Ch18  Ch19  Ch20  Ch21  Ch22  Ch23  Ch24  Ch25  Ch26  Ch27  Ch28  Ch29  Ch30  Ch31  Ch32  Ch33  Ch34  Ch35  Ch36  Ch37  Ch38  Ch39
%  1   2  3      4      5       6           7        8     9    10   11  12   13 
%  


clear all
close all

home_dir = pwd;
addpath(home_dir)
%% INPUT FILES
cd 'C:\LC\Private\dominiquef\Projects\4264_Vedanta_TEM\Processing\HTS_data'
timefile='HTS_channel_times.chn';

% Data components
datafile = 'HTS_data_loop3_p125Hz_Zcomp.xyz';
datafile_Y = [];
datafile_X = [];

loopfile='Loop3_HTS_densify_drape.tx';

current_dir = 'counterclockwise';

%% Specifications
isdataLTS=0; %LTS SQUID data or HTS SQUID data?

Lasttime    = 45; %number of channels to drop from end
Firstime   = 13;

ChanSpacing=2;

% First time chanel column in file
index = 8;

% Uncertainties
ErrorPerc=0.05; %percent

ErrorFloorHTS=0.5; %minimum erro floor allowed
ErrorFloorLTS=0.25;%pT/A



%units of B are pT/A
%convert by multiplying by 1e-12 to T and then divide by mu_0 to get A/m
mu_0=4*pi*1e-7;


%% Load Times
times = load(timefile);
%convert to seconds
times=times*1e-3; %channel times
TimingMark=0.5*1e-3;% ramp shutoff at 0.5 ms  

% times=times+TimingMark;



%% Load Data File
data=load(datafile);
ndv = data(:,7)==-1e+33;
data = data(ndv==0,:);

% if isempty(datafile_Y)==0
%     Dy=load(datafile_Y);
%     ndv = Dy(:,7)==-1e+33;
%     Dy = Dy(ndv==0,:);
% end
fidout=fopen([datafile,'_ALL_UBC_format.out'],'w');


% D = load(datafile_Y);
% ndv = D(:,8)==-1e+33;
% D = D(ndv==0,:); 
% 
% Dz =[D(D(:,3)==3,1:2) D(D(:,3)==3,4:end)];
% Dy =[D(D(:,3)==1,1:2) D(D(:,3)==1,4:end)];



%% Load Loop
L_in = load(loopfile);

% Check loop
% 'flag': 'clockwise' or 'counterclockwise'
L = checkloop(L_in,current_dir);

%Specify earliest time gate to use



% for ii=1:1:length(times)
% if times(ii,1)*1000 > MinTime
%     break;
% end
% end
% 
% DropEarly=ii; %number of early channels to drop

%number to increment the channel spacing i.e. 2 means uese every other
%channel


%% Set ERROR floor

if isdataLTS==1
    
    fprintf('Low Temp SQUID data error floor used of %f\n',ErrorFloorLTS);
    ErrorFloor=ErrorFloorLTS * 1e-12 /mu_0;
    
else
 
    ErrorFloor=ErrorFloorHTS * 1e-12 /mu_0;   

end


%% Convert B-field to H-field


data(:,index:end)=(data(:,index:end)*1e-12)/mu_0;
% Dy(:,index:end)=(-Dy(:,index:end)*1e-12)/mu_0;

% Pre-allocate memory
ErrorFloor_z = zeros(size(data,2),1);
% ErrorFloor_y = zeros(size(Dy,2),1);

% Assign error floor
for ii = index : size(data,2)
    
    ErrorFloor_z(ii) = max( ErrorPerc * median(abs(data(:,ii))) , ErrorFloor );
    
%     ErrorFloor_y(ii) = max( 2 * ErrorPerc * median(abs(Dy(:,ii))) , ErrorFloor );
    
end

%% Write data file header

ig=99;


% header
fprintf(fidout,'IGNORE %i\n',ig);
fprintf(fidout,'\n');

fprintf(fidout,'N_TRX 1\n');
fprintf(fidout,'\n');


%% write transmitter
fprintf(fidout,'TRX_ORIG\n');
S=size(L);
fprintf(fidout,'%i\n',S(1));
for ii=1:1:S(1)
    fprintf(fidout,'%6.2f %6.2f %6.2f\n',L(ii,:));
end

fprintf(fidout,'\n');


%% write data    
ST=size(times);

%% Check for duplicates
% if isempty(datafile_Z) == 0
ii = 1;
while ii <  size(data,1)
    
    xmatch = data(ii,1) == data(:,1);
    ymatch = data(ii,2) == data(:,2);
%     zmatch = D(ii,4) == D(:,4);
    
    temp = (xmatch .* ymatch) == 1;
    [dupl] = find(temp);
    
    if length(dupl)>1
        fprintf('Found duplicate station location. Data will be average\n')
    end
    % Average the values with same location, then remove replicate lines
    data(dupl(1),:) = mean(data(temp,:),1);
    
    % Collapse replicate lines
    temp(dupl(1))=0;
    data = data(temp==0,:);
    ii = ii + 1;
    
end
save([datafile '.avg'],'-ascii','data');
% end

% if isempty(datafile_Y) == 0
% % Check for duplicates in Y component
% ii = 1;
% while ii <  size(Dy,1)
%     
%     xmatch = Dy(ii,1) == Dy(:,1);
%     ymatch = Dy(ii,2) == Dy(:,2);
% %     zmatch = D(ii,4) == D(:,4);
%     
%     temp = (xmatch .* ymatch) == 1;
%     [dupl] = find(temp);
%     
%     if length(dupl)>1
%         fprintf('Found duplicate station location. Data will be average\n')
%     end
%     % Average the values with same location, then remove replicate lines
%     Dy(dupl(1),:) = mean(Dy(temp,:),1);
%     
%     % Collapse replicate lines
%     temp(dupl(1))=0;
%     Dy = Dy(temp==0,:);
%     ii = ii + 1;
%     
% end
% save([datafile_Y '.avg'],'-ascii','Dy');
% 
% end

%%
S=size(data);

numchans_indatafile=S(2)-(index-1);
numchans_intimesfile=ST(1);


if numchans_indatafile ~= numchans_intimesfile
     fprintf('number of times in data file and number of times in channel file do not match\n');

     keyboard;
     
else
    numchans=numchans_indatafile;
end

jj=Firstime:ChanSpacing:Lasttime;  
n_time=length(jj);
n_recv=S(1);

%data header

fprintf(fidout,'N_RECV %i\n',n_recv);
fprintf(fidout,'N_TIME %i\n',n_time);
    fprintf(fidout,'\n');

    %data records

for ii=1:1:S(1)       
            
    for jj=Firstime:ChanSpacing:Lasttime                      
    fprintf(fidout,'%7.2f  %7.2f %7.2f ',data(ii,1),data(ii,2),data(ii,3));
    fprintf(fidout,'%7e ',times(jj,1));   
    
    % Look for coordinate match between component
%     xmatch = Dz(ii,1) == Dy(:,1);
%     ymatch = Dz(ii,2) == Dy(:,2);
%     zmatch = D(ii,4) == D(:,4);
    
%     temp = (xmatch .* ymatch) == 1;
    
    for mm=1:1:10
        fprintf(fidout,'%2.0f ',ig);
    end
    
%     if sum(temp) == 1
%         fprintf(fidout,'%e %e ',Dy(temp,jj+(index-1)),abs(Dy(temp,jj+(index-1))*2*ErrorPerc)+ErrorFloor);
%     
%     else
%         fprintf(fidout,'%2.0f ',ig);    
%     end
    
    fprintf(fidout,'%e %e ',data(ii,jj+(index-1)),abs(data(ii,jj+(index-1))*ErrorPerc)+ErrorFloor);

        for mm=1:1:6
        fprintf(fidout,'%2.0f ',ig);
        end

    fprintf(fidout,'\n');
    end
end


fclose all;