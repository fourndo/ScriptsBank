%% E3D_Total_Field_conv
% Writte for DIGHEM data measured in ppm
%
% Load a PRIMARY.pre file (FREESPACE), then reads in aE3D.obs (ppm) file 
% and write a new obsfile with data converted to Total Field 
% In-phase      = ppm * 1e-6 * Hp + Hp
% Quadrature    = ppm * 1e-6 * Hp
%
% Hp is computed using E3Dfwd and a freespace model (1e-8)
% 
% Writte for DIGHEM data measured in ppm
%
% Written by: D.Fournier
% Last Update: September 17th, 2014

clear all
close all

% addpath C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB

%% USER INPUT
work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\3D\DO27\Inv40_newUncert';
obsfile     = 'DIGHEMdwn.txt';
predfile    = 'e3d_data_FREESPACE.txt';
outfile     = 'Total_field_data.dat';



%% SCRIPT STARTS HERE

% Start GIFTools Project
EMproj = GIFproject();

% Load data
obs = E3Dinversion.readDobs(EMproj,[work_dir '\' obsfile]); 

% Extract transmitters
tx = obs.getTx;
rx = obs.getData(:,1:5);

% Load predicted 
% !! need to run e3dfwd with freespace model first !!
E3Dinversion.readPred(EMproj,obs,[work_dir '\' predfile]);
pred = EMproj.getItem(2); 
pred.setName('data_FREESPACE'); 
pred.setioHeadInd([3   4   5   6   8  10  12  14  16  18  20  22  24  26  28   1   2;],{'X','Y','Z','ExR','ExI','EyR','EyI','EzR','EzI','HxR','HxI','HyR','HyI','HzR','HzI','TXID','FREQUENCY'}); 

freq = obs.getData(:,'FREQUENCY');
%% Re lines in obsfile
% Assign Uncertainties by frequency and real/imag
% Might need to adapt if uncertainties already exist
uncert_R = zeros( length(freq),1 );
uncert_I = zeros( length(freq),1 );

uncert_R( freq==900 ) = 1;%(std(data(data(:,1)==900,end-3))) * uncert_pct;
uncert_I( freq==900 ) = 1;%(std(data(data(:,1)==900,end-1))) * uncert_pct;

uncert_R( freq==7200 ) = 2;%(std(data(data(:,1)==7200,end-3))) * uncert_pct;
uncert_I( freq==7200 ) = 2;%(std(data(data(:,1)==7200,end-1))) * uncert_pct;

uncert_R( freq==56000 ) = 3;%(std(data(data(:,1)==56000,end-3))) * uncert_pct;
uncert_I( freq==56000 ) = 3;%(std(data(data(:,1)==56000,end-1))) * uncert_pct;

%% Convert from ppm to Total Field Hz only for
% Convert data and uncertainties to total field
Hz_R = obs.getData(:,'HzR');
Hz_I = obs.getData(:,'HzI');

% Get primary field
Ho_R = pred.getData(:,'HzR');

% Convert all to primary
Hz_R = Hz_R .* Ho_R * 1e-6; 
Hz_I = Hz_I .* Ho_R * 1e-6;

% Add primary to in-phase data only
Hz_R = Hz_R + Ho_R;

% Convert uncertainties
uncert_R = uncert_R .* Ho_R * 1e-6; 
uncert_I = uncert_I .* Ho_R * 1e-6;


%% Write data to file
fid = fopen([work_dir '\' outfile],'w');

fprintf(fid,'! Export from E3D_Total_Field_conv: D.Fournier\n');
fprintf(fid,'IGNORE NaN\n\n');
fprintf(fid,'N_TRX %i\n', size(tx,1) );


for ii = 1 : size(rx,1);
    
    

    tx_spec = tx{rx(ii,1)}.getData;

    fprintf(fid,'TRX_LOOP\n');
    fprintf(fid,'%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n',tx_spec(1:6));
    fprintf(fid,'FREQUENCY %12.8e\n',freq(ii));
    fprintf(fid,'N_RECV %i\n', 1 );


    % Write receivers   
    fprintf(fid,'%12.8e %12.8e %12.8e ',rx(ii,3:5));

    for ll = 1 : 22

        fprintf(fid,'NaN ');

    end

    fprintf(fid,'%12.8e %12.8e %12.8e %12.8e',...
        Hz_R(ii),uncert_R(ii), Hz_I(ii),uncert_I(ii));

    fprintf(fid,'\n \n');
            

    
end

fclose(fid);
