% Script to assign uncertainties and format ZTEM data
clear all
close all

addpath '..\FUNC_LIB'
%% INPUT PARAMETERS
work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\ZTEM\FLIN';
datafile = 'FLIN_ZTEM_ROT_full.dat';

dsep = '\';

% Define column index for the following components
freq = [30 45 90 180 360 720];

% Base station location
xyz0 = [311750 6075530 340];
%% SCRIPT STARTS HERE
d=load([work_dir dsep datafile]);

xyz     = d(:,1:3);
data    = d(:,4:end); 

ndv = -99999;
%% OPTIONAL FILTER DATA
radius = 50;

figure;
scatter(xyz(:,1),xyz(:,2),'b*'); hold on

[indx] = Filter_xy(xyz(:,1),xyz(:,2),radius);

data = data(indx==1,:);
xyz = xyz(indx==1,:);

scatter(xyz(:,1),xyz(:,2),'ro')

ndata = size(data,1);

% 1  X
% 2  Y
% 3  Z
% 4  xip_30Hz
% 5  xip_45Hz
% 6 xip_90Hz
% 7 xip_180Hz
% 8 xip_360Hz
% 9 xip_720Hz
% 10 xqd_30Hz
% 11 xqd_45Hz
% 12 xqd_90Hz
% 13 xqd_180Hz
% 14 xqd_360Hz
% 15 xqd_720Hz
% 16 yip_30Hz
% 17 yip_45Hz
% 18 yip_90Hz
% 19 yip_180Hz
% 20 yip_360Hz
% 21 yip_720Hz
% 22 yqd_30Hz
% 23 yqd_45Hz
% 24 yqd_90Hz
% 25 yqd_180Hz
% 26 yqd_360Hz
% 27 yqd_720Hz


%% assign standard deviations
uncert = zeros(ndata,size(data,2));

for ii = 1 : size(data,2)
    uncert(:,ii) = 1e-1 + 0.05*abs(data(:,ii));

%     uncert(:,ii) = 0.5*repmat(min(abs(data(data(:,ii)~=ndv,ii))),ndata,1) + 0.02*abs(data(:,ii));
    uncert(data(:,ii)==ndv,ii) = ndv;
end

%% Optional cropping

% d_crop=[];
% for i=1:s_d(1)
%     if d(i,1)>406000
%         d_crop=[d_crop; d(i,:)]; %#ok<AGROW>
%     end
%     
% end
% 
% d=d_crop;
% 
% s_d=size(d)

%% Geotech 2 UBC sign Convension

% UBC_xip=-yip;
% UBC_xqd=-yqd;
% UBC_yip=-xip;
% UBC_yqd=-xqd;


% T_30 =[d(:,1:3) -d(:,20) d(:,44)  -d(:,26) d(:,50)  -d(:,8) d(:,32) -d(:,14) d(:,38)];
% T_45 =[d(:,1:3) -d(:,21) d(:,45)  -d(:,27) d(:,51)  -d(:,9) d(:,33) -d(:,15) d(:,39)];
% T_90 =[d(:,1:3) -d(:,22) d(:,46)  -d(:,28) d(:,52)  -d(:,10) d(:,34) -d(:,16) d(:,40)];
% T_180=[d(:,1:3) -d(:,23) d(:,47)  -d(:,29) d(:,53)  -d(:,11) d(:,35) -d(:,17) d(:,41)];
% T_360=[d(:,1:3) -d(:,24) d(:,48)  -d(:,30) d(:,54)  -d(:,12) d(:,36) -d(:,18) d(:,42)];
% T_720=[d(:,1:3) -d(:,25) d(:,49)  -d(:,31) d(:,55)  -d(:,13) d(:,37) -d(:,19) d(:,43)];
% 
% 
% save T_30.dat T_30 -ascii
% save T_45.dat T_45 -ascii
% save T_90.dat T_90 -ascii
% save T_180.dat T_180 -ascii
% save T_360.dat T_360 -ascii
% save T_360_2.dat T_720 -ascii

%% Write obs file
fid = fopen([work_dir dsep 'ZTEM.obs'],'w');

fprintf(fid,'DATATYPE MTT\n');
fprintf(fid,'IGNORE %i\n',ndv);

for ii = 1 : length(freq)
    
    fprintf(fid,'\n');
    fprintf(fid,'%8.5e ! Frequency \n',freq(ii));
    fprintf(fid,'%i\n',ndata +1 );
    fprintf(fid,'%8.5e %8.5e %8.5e %i %i %i %i %i %i %i %i\n',xyz0(1),xyz0(2),xyz0(3),ndv,ndv,ndv,ndv,ndv,ndv,ndv,ndv);
    
    for jj = 1 : ndata
        
        fprintf(fid,'%8.5e %8.5e %8.5e ',xyz(jj,:));
        
        %Print real X
        if data(jj,1+(ii-1)) == ndv
            fprintf(fid,'%i %i ',ndv,ndv);
            
        else
        
            fprintf(fid,'%8.5e %8.5e ',-data(jj,1+(ii-1)),uncert(jj,1+(ii-1)));
        end
        
        %Print imag X
        if data(jj,7+(ii-1)) == ndv
            fprintf(fid,'%i %i ',ndv,ndv);
            
        else
            fprintf(fid,'%8.5e %8.5e ',-data(jj,7+(ii-1)),uncert(jj,7+(ii-1)));
        end
        
        %Print real Y
        if data(jj,13+(ii-1)) == ndv
            fprintf(fid,'%i %i ',ndv,ndv);
            
        else
        
            fprintf(fid,'%8.5e %8.5e ',-data(jj,13+(ii-1)),uncert(jj,13+(ii-1)));
        end
        
        %Print imag Y
        if data(jj,19+(ii-1)) == ndv
            fprintf(fid,'%i %i ',ndv,ndv);
            
        else 
        
            fprintf(fid,'%8.5e %8.5e',-data(jj,19+(ii-1)),uncert(jj,19+(ii-1)));
        
        end
        fprintf(fid,'\n');
        
    end
    
end

fclose(fid);