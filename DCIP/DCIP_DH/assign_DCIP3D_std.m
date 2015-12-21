% function assign_DCIP3D_std(work_dir,obsfile,pct,floor)
% Function assign_DCIP3D_std(obsfile,pct,floor)
% Load a UBC-DCIP3D observation file and assign error to each datum
% such as: error = |datum| * pct/100 + floor
% 
% INPUT
% work_dir: working directory for the obsfile and output result
% obsfile : observation file (DCIP3D format)
% pct     : percentage of each datum assign as error (pct/100)
% floor   : minimum error floor
%
% Written by: Dom Fourner
% Last update: 2013-12-08

%% FOR DEV ONLY
close all
clear all

obsfile = 'dc3d.dat';
work_dir = 'C:\Projects\4180_Wallbridge_Wisner_IP\Data\DataProcessing\DHIP\FWR\Inv_synthetic';
pct = 5;
floor = 1e-4;

%/////////////////////%
%% SCRIPT START HERE %%
%/////////////////////%

% Load obsfile
fid = fopen([work_dir '\' obsfile],'r');
wrt2file = fopen([work_dir '\DCIP3D_' obsfile(1:end-4) num2str(pct) 'pct_' num2str(floor) 'floor.dat'],'w');

% Cycle through all the lines and add data to outfile
line=fgets(fid);

while line~=-1
    
    data = str2num(line);
    
    % Find transmitter line
    if isempty(data)== 0

        % Get the number of receiver data
        nrx = data(end);
        fprintf(wrt2file,'%s',line);
        
        % Cycle down the rx and assign error
        for jj = 1:nrx
            
            line=fgets(fid);
            data = str2num(line);
            
            std = abs(data(end)) * pct/100 + floor;
            
            for jj = 1 : length(data)
                
                fprintf(wrt2file,'%12.8e ',data(jj));
                
            end
            
            fprintf(wrt2file,'%12.8e\n',std);
            
        end
        
    else
        
        fprintf(wrt2file,'%s\n',line);
        
    end
    
    line=fgets(fid);
    
end
    
    
        
            