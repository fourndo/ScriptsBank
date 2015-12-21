clear all
close all

% Load DIGHEM data and assign uncertainties
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\GIF_project';

% Project name
proj = 'TKC_DIGHEM_Coarse';
cd(work_dir);

load(proj)
project = GIFtools(obj);

%% Load data from ASCII
% FEMdata.importCSV(project,[project.getWorkDir,'\.\TKC_DIGHEM_edt.txt']); 
% item = project.getItem(1); item.sortImportedData({'X_NAD27','Y_NAD27','Bird'},[900;7200;56000;],{'CPQ900','CPI900';'CPQ7200','CPI7200';'CPQ56k','CPI56k'},{'HzI';'HzR'}); 
% item = project.getItem(1); item.setName('dighem'); 
% item = project.getItem(1); item.setioHeadInd([0;1;2;3;4;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;6;0;5;0;]); 

%% Assign uncertainty per frequency
% item = project.getItem(1); item.addDataTypeStn({'HzI'},{'HzIuncert'},[1;2;4;],[0.1;0.1;0.1;]); 
% item = project.getItem(1); item.addDataTypeStn({'HzR'},{'HzRuncert'},[1;2;4;],[0.1;0.1;0.1;]); 
% item = project.getItem(1); item.setioHeadInd([0;1;2;3;4;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;6;8;5;7;]); 

%% grid Topo and get topo at our measurement points
% 
% % Load topo
% gridTopo = TriScatteredInterp([NAD27x NAD27y], Elevatn);
% Topo     = gridTopo([item.getData(:,'X_NAD27') item.getData(:,'Y_NAD27')]);
% 
% % add topo and Bird elevation columns to our data
% item.addColumn('Topo',Topo); % topography at measurement points
% item.addColumn('Bird_Elevatn',Topo+item.getData(:,'Bird')); % topo + bird (bird position in elevation)


%item = project.getItem(1);
% data = item.getData;
% 
% item.replaceColumn(4,data(:,4) + 423);

%% Assign bearing direction for trasnmitter offset
% item = project.getItem(1);
% item.calculateBearing;

%% Create transmitters
% item.addCoplanarDIGHEMsources([],[],'bearing')

%% Flip the sign of Imaginary
% item = project.getItem(1);
% d = item.getData;
% 
% item.replaceColumn(5,-d(:,5));

%% Change uncertainties
item = project.getItem(21);
data = item.getData;
item.copy;

HzR_uncrt = zeros(size(data,1),1);
HzI_uncrt = zeros(size(data,1),1);

index = data(:,2)==900;

HzR_uncrt(index) = abs(data(index,6) ) * 0.0 + 1;%max([0 0.05*abs(min(data(index,6)))]);
HzI_uncrt(index) = abs(data(index,8 )) * 0.0 + 1;%max([0 0.05*abs(max(data(index,8)))]);

index = data(:,2)==7200;

HzR_uncrt(index) = abs(data(index,6) ) * 0.0 + 2;%max([0 0.1*abs(min(data(index,6)))]);
HzI_uncrt(index) = abs(data(index,8) ) * 0.0 + 2;%max([0 0.1*abs(max(data(index,8)))]);

index = data(:,2)==56000;

HzR_uncrt(index) = abs(data(index,6) ) * 0.0 + 4;%max([0 0.1*abs(min(data(index,6)))]);
HzI_uncrt(index) = abs(data(index,8) ) * 0.0 + 4;%max([0 0.1*abs(max(data(index,8)))]);

item = project.getItem(end);
item.replaceColumn([7 9],[HzR_uncrt HzI_uncrt]);

%% Load pred data and modify to ppm
item_ppm = project.getItem(22);
Ho_R = item_ppm.getData(:,16);

item = project.getItem(end);
data = item.getData;

% Find fields smaller and equal to 0 and replace by nan
% index = data(:,2)==900;
% 
% data(index,6) = data(index,6)+2;
% 
% index = data(:,2)==7200;
% 
% data(index,6) = data(index,6)+1;

% data(index,9) = nan;

% Find fields smaller and equal to 0 and replace by nan
index = data(:,8)>-1;
% 
% data(index,6) = 0.5;
data(index,8) = -0.01;

% Find fields smaller and equal to 0 and replace by nan
index = data(:,6)<=0;
% 
% data(index,6) = 0.5;
data(index,6) = 0.01;

% Convert all to primary
data(:,6) = data(:,6) .* Ho_R * 1e-6; 
data(:,7) = data(:,7) .* Ho_R * 1e-6; 
data(:,8) = data(:,8) .* Ho_R * 1e-6; 
data(:,9) = data(:,9) .* Ho_R * 1e-6; 

% Add primary to in-phase data only
data(:,6) = data(:,6) + Ho_R;

item = project.getItem(end);
item.replaceColumn(6:9,data(:,6:9));

%% After loading predicted, compute convert to ppm and compute residual
item = project.getItem(19);
Ho_R = item.getData(:,16);

item = project.getItem(18);
obs = item.getData;

item = project.getItem(end);
pred = item.getData;

% Hz_R = (pred(:,4) - Ho_R);
Hz_R = (pred(:,16) - Ho_R) ./ Ho_R * 1e+6;
Hz_I =  pred(:,17) ./ Ho_R * 1e+6;

norm_res_Hz_R = (obs(:,6)  - Hz_R )./obs(:,7);
norm_res_Hz_I = (obs(:,8) -  Hz_I )./obs(:,9);

% item.addColumn('Hz_R_ppm',Hz_R); 
% item.addColumn('ObsHz_R_ppm',obs(:,6));
% item.addColumn('Hz_I_ppm',Hz_R); 
% item.addColumn('ObsHz_I_ppm',obs(:,8)); 

item.addColumn('Hz_R',Hz_R); 
item.addColumn('ObsHz_R',obs(:,6));
item.addColumn('Hz_I',Hz_I); 
item.addColumn('ObsHz_I',obs(:,8)); 

item.addColumn('norm_res_ObsHz_R',norm_res_Hz_R); 
item.addColumn('norm_res_ObsHz_I',norm_res_Hz_I); 


%% Write to XYZ format

out_mat = zeros(size(item.tx,1),21);

d_pred = item.getData;

indx = d_pred(:,2)==900;
out_mat(:,1) = d_pred(indx,3);   % X-coord.
out_mat(:,2) = d_pred(indx,4);   % Y-coord.
out_mat(:,3) = d_pred(indx,5);   % Z-coord.

indx = d_pred(:,2)==900;
% out_mat(:,4) = d_pred(indx,19);   % Obs-900-R
out_mat(:,5) = d_pred(indx,18);   % Pre-900-R
% out_mat(:,6) = d_pred(indx,22);   % Res-900-R
% out_mat(:,7) = d_pred(indx,21);   % Obs-900-R
out_mat(:,8) = d_pred(indx,20);   % Pre-900-R
% out_mat(:,9) = d_pred(indx,23);   % Res-900-R

indx = d_pred(:,2)==7200;
% out_mat(:,10) = d_pred(indx,19);   % Obs-7200-R
out_mat(:,11) = d_pred(indx,18);   % Pre-7200-R
% out_mat(:,12) = d_pred(indx,22);   % Res-7200-R
% out_mat(:,13) = d_pred(indx,21);   % Obs-7200-R
out_mat(:,14) = d_pred(indx,20);   % Pre-7200-R
% out_mat(:,15) = d_pred(indx,23);   % Res-7200-R

indx = d_pred(:,2)==56000;
% out_mat(:,16) = d_pred(indx,19);   % Obs-56k-R
out_mat(:,17) = d_pred(indx,18);   % Pre-56k-R
% out_mat(:,18) = d_pred(indx,22);   % Res-56k-R
% out_mat(:,19) = d_pred(indx,21);   % Obs-56k-R
out_mat(:,20) = d_pred(indx,20);   % Pre-56k-R
% out_mat(:,21) = d_pred(indx,23);   % Res-56k-R

save([work_dir '\Obs_vs_pred_RENAME.dat'],'-ascii','out_mat')

%% Load octree code and create cell face to compensate for alpha
[G] = mesh.getCenterGradientMatrix;
aa = ones(size(G,2),1);
aa(1:2:end) = 0;

gradaa = G*aa;

gradaa(gradaa~=0) = 1./gradaa(gradaa~=0); 
gradaa(gradaa==0) = 1;
gradaa = gradaa.^2;
save([work_dir '\FweightV.dat'],'-ascii','gradaa');
