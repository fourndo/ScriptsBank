% Data Tiled Merge
%
% INPUT: 
% work_dir: Working directory
% OPTIONAL
% meshfile : Final mesh
% 
%
% OUTPUT:
% Final merged model
% Author: D.Fournier
%
% Last update: June 9th, 2015
%

clear all
close all


%% INPUTS PARAMETERS
% Load files
work_dir ='C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\CMI_Tiles\ALL_Tiles';

% Observed file
obs_file = '..\..\Obs_Paul_Lake_SUB_5pc_5nT_DETREND.dat';

% Tile file
tile_file = '..\..\Tiles_50m_Tight.dat';

% Maximum search radius (m), recommended x2 smallest cell size
rangemax = 500;

%No data value (usually -1(Mag), -100(Grav), -99999(Gocad), 1e-008(DC))
ndv=-100;

dsep = '\';
%% Generate list of input files

filein = ls(work_dir);
filenum = [];
for ii = 1 : size(filein,1)
    
    temp1 = regexp(filein(ii,:),'\d*','match');
    temp2 = regexp(filein(ii,:),'[_]','split');
    if ~isempty(temp1) && size(temp2,2) > 1
        
        if strcmp(strtrim(temp2{2}),'MVI.pre')
        
            filenum = [filenum str2num(temp1{:})];
        
        end
        
    end
end

tt = unique(filenum);
ntiles = length(tt);

%% Load in obs

tiles = load([work_dir '\' tile_file]);

% ntiles = size(tiles,1);

% Input files
for ii = 1 : ntiles
    
%     if ~isempty(regexp(mod_list(ii,:),'pre','match'))
%         count = count + 1;
%         
%         in_file = [strtrim(mod_list(ii,:));
        
        pre_file{ii}=[work_dir dsep 'Tile' num2str(tt(ii)) '_MVI.pre'];

        
%     end
    
end

%% \\\\\\\\\\\\\\\\\ %%
% SCRIPT STARTS HERE  %
%%\\\\\\\\\\\\\\\\\\ %%

% Load observed data
[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd] = read_MAG3D_obs([work_dir dsep obs_file]);

% Round location to avoid rounding error
% Obsx = round(Obsx);
% Obsy = round(Obsy);
% Obsz = round(Obsz);

% Load predicted files for all the tiles
for ii = 1 : ntiles
        
    [~, ~, ~, ~, obsx{ii}, obsy{ii}, ~, d{ii}, ~] = read_MAG3D_obs(pre_file{ii});
    
    cntr_x = mean(obsx{ii});
    cntr_y = mean(obsy{ii});
    
    % Create distance weights
    w{ii} = ( ( cntr_x - obsx{ii} ).^2 + ( cntr_y - obsy{ii} ).^2 ).^(-0.5);
    
end

%% Create final predicted data
pred = zeros(size(data));
denom = zeros(size(data));

count = 1;
progress = -1;
tic
for ii = 1 : ntiles

    indx = Obsx >= tiles(tt(ii),1) & Obsx <= tiles(tt(ii),3) & Obsy >= tiles(tt(ii),2) & Obsy <= tiles(tt(ii),4);
    
    pred(indx) = pred(indx) + d{ii}.*w{ii};
    denom(indx) = denom(indx) + w{ii};

    

    % Monitor progress and print to screen
     d_iter = floor(ii/ntiles*20);
    if  d_iter > progress

        fprintf('Merged %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            

        
end

pred(denom>0) = pred(denom>0)./denom(denom>0);

plot_TMI(Obsx,Obsy,data,pred,wd,'Paul Lake')

% Only output non-zero pred
indx = pred~=0;
%% Save predicted to file
write_MAG3D_TMI([work_dir dsep 'Merged_MVI_TMI.pre'],H,I,Dazm,Obsx(indx),Obsy(indx),Obsz(indx),pred(indx),wd(indx))
    