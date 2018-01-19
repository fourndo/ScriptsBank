%% E3D_Total_Field_conv
% Load a PRIMARY.pre file (FREESPACE), then reads in aE3D.obs (ppm) file 
% and write a new obsfile with data converted to Total Field 
% In-phase      = ppm * 1e-6 * Hp + Hp
% Quadrature    = ppm * 1e-6 * Hp
%
% Hp is computed using E3Dfwd and a freespace model (1e-8)
% 
% Written by: D.Fournier
% Last Update: September 17th, 2014

clear all
close all

addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

%% USER INPUT
work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\Real_topo\Inv12_PPMx2';
obsfile     = 'dighem_dwns_20m_ALLfreq_ppm_Ztopo.dat';
predfile    = 'FWR_dighem_dwns_20m_ALLfreq_TField_Ztopo.dat';
outfile     = 'Data_TotalField_0pct_varflr_PPMx2.dat';

pct_e = 0.0;
%% SCRIPT STARTS HERE
% Load data
[trx,data] = load_E3D_obs([work_dir '\' obsfile]);

% Scale data by 1/2 for co-planar_x2
data(:,[25 27]) = data(:,[25 27])*2.0;

% Load pred
[~,pred] = load_E3D_obs([work_dir '\' predfile]);
Hp = pred(:,1:16);

% Hp = load_E3D_pred([ work_dir '\' predfile]);


%% Write obsfile from FWR modelled data
% freq = unique(pred(:,1));
% 
% % Add noise
% noise = 0.001;
% data = pred(:,end) +  pred(:,end) .* noise .* randn(size(pred,1),1);
% uncert = abs(pred(:,end)) * noise + noise * std((pred(:,end)));
% 
% write_e3d_obs([work_dir '\Synthetic_small.dat'],freq,trx,pred(:,1:3),[data uncert]);

%% Re lines in obsfile
filein = fopen([work_dir '\' obsfile],'r');
fileout = fopen([work_dir '\' outfile],'w');

freq = unique(data(:,1));

% Assign Uncertainties by frequency

%uncert_pct = 0.05;
% obs(:,end-1:end) = obs(:,end-1:end) +  obs(:,end-1:end) * 1e-6 .* 0.1 .* randn(size(obs,1),2);
uncert(1,1) = 1;%(std(data(data(:,1)==900,end-3))) * uncert_pct;
uncert(1,2) = 1;%(std(data(data(:,1)==900,end-1))) * uncert_pct;

uncert(2,1) = 2;%(std(data(data(:,1)==7200,end-3))) * uncert_pct;
uncert(2,2) = 2;%(std(data(data(:,1)==7200,end-1))) * uncert_pct;

uncert(3,1) = 3;%(std(data(data(:,1)==56000,end-3))) * uncert_pct;
uncert(3,2) = 3;%(std(data(data(:,1)==56000,end-1))) * uncert_pct;

% uncert = [uncert abs(std(obs(1:2,end))) * 1e-1;]; %+ 0.01 * std((obs(:,end)));

indx = ones(size(data,1),1);

for ii = 1 : size(data,1)
    
    if data(ii,25) <=0
        
        data(ii,25:26) = NaN;
        
    else
        
        data(ii,26) = abs(data(ii,25))*pct_e + uncert( freq==data(ii,1) , 1 );
    
    end
    
    if data(ii,27) >=0
        
        data(ii,27:28) = NaN;
        
    else
        
        data(ii,28) = abs(data(ii,27))*pct_e + uncert( freq==data(ii,1) , 2 );
    
    end
    
    if isnan(data(ii,[25 27]))
        
        indx(ii) = 0;
        
    end
    
    % Convert data and uncertianties to total field
    data(ii,25:28) = data(ii,25:28) * 1e-6 * Hp(ii,end-1);

    % Add primary to in-phase data only
    data(ii,25) = data(ii,25) + Hp(ii,end-1);
    
end

trx = trx(indx==1,:);
tid = kron(1:size(trx,1),ones(1,3));

data = data(indx==1,:);
tid = tid(indx==1);
write_e3d_obs([work_dir '\' outfile],trx,data)


% line = fgets(filein);
% count = 1;
% while line~=-1
%     
%     if isempty(regexp(line,'N_RECV','match'))==0
%         
%         fprintf(fileout,'%s',line);
%         nrecv =regexp(line,'\s\d*','match');
%         nrecv = str2num(nrecv{1});
%         count_recv = 0;
%         while count_recv < nrecv
%             line = fgets(filein);
%             
%             if isempty(str2num(line))==0
%                 
%                 temp = str2num(line);
%                 
% %                 if temp(end-3) == 0 || temp(end-1)==0
% %                     
% %                     temp(end-3:2:end-1) = NaN;
% %                     temp(end-2:2:end) = NaN;
% %                     
% %                 else
% 
%                     % Add floor
%                     temp(end-2) = abs(temp(end-3))*pct_e + uncert( freq==data(count,1) , 1 );%uncert;
%                     temp(end)   = abs(temp(end-1))*pct_e + uncert( freq==data(count,1) , 2 );%uncert;
%                     
%                     % Convert data and uncertianties to total field
%                     temp(end-3:end) = temp(end-3:end) *1e-6 * Hp(count,end-1);
%                     
%                     % Add primary to in-phase data only
%                     temp(end-3) = temp(end-3) + Hp(count,end-1);
%                     
% %                 end
%                 
%                 for kk = 1 : length(temp); 
%                     fprintf(fileout,'%12.8e ',temp(kk));
%                 end
%                 
%                 fprintf(fileout,'\n');
%                 
%                 count = count + 1;
%                 count_recv = count_recv+1;
%             end
%             
%         end
%         
%     else
%         
%         fprintf(fileout,'%s\n',line);
%         
%     end
%     
%     line = fgets(filein);
%     
% end
%     
%       
% fclose all;