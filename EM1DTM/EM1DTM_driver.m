% FOR DEV ONLY
% Temporary function to format data to GIFTool format
% MUST BE REFORMATED FOR EVERY FILE

clear all
close all

addpath '.\Functions';

work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\VTEM\1D';

dsep = '\';

%% Load input file
[meshfile,obsfile,topofile,nullfile,m_con,con_ref,m_sus,sus_ref,alpha_con,alpha_sus,beta,cooling,target,bounds,mtype,interp_n,interp_r,interp_s] = EM1DTM_read_inp([work_dir dsep 'EM1DTM_LC.inp']);

%% Load models and parameters
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);

[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

mcell = (length(zn)-1)*(length(xn)-1)*(length(yn)-1);


% Create or load reference cond model
if ischar(m_con)==1
    
    m_con = load([work_dir dsep m_con]);
    HSflag = 1;
elseif isempty(m_con)

    % If m_con is empty, then start with HS
    m_con = ones(mcell,1)*1e-8;
    HSflag = 0;
    
else
    
    m_con = ones(mcell,1)*m_con;
    HSflag = 1;
end

% Create or load starting cond model
if ischar(con_ref)==1
    
    con_ref = load([work_dir dsep con_ref]);
    
else
    
    con_ref = ones(mcell,1)*con_ref;
    
end

% Create or load reference cond model
if ischar(m_sus)==1
    
    m_sus = load([work_dir dsep m_sus]);
    
else
    
    m_sus = ones(mcell,1)*m_sus;
    
end

% Create or load starting cond model
if ischar(sus_ref)==1
    
    sus_ref = load([work_dir dsep sus_ref]);
    
else
    
    sus_ref = ones(mcell,1)*sus_ref;
    
end


%% Load topography
if ~isempty(topofile) && isempty(nullfile)
    
    topo = read_UBC_topo([work_dir dsep topofile]);
    [nullcell,temp,temp] = topocheck(xn,yn,zn,topo);
    save([work_dir dsep 'nullcell.dat'],'-ascii','nullcell');

elseif isempty(topofile) && ~isempty(nullfile)
    
    nullcell = load([work_dir dsep nullfile]);
else
    
    nullcell = ones(mcell,1);
    
end
%% Reformat data in array structure
% Last entry specifies the minimum distance between points
% sort_EM1DTM will downsample the data using a minimum distance
% as specified by the last input

%% VTEM SURVEY
% NEED TO WRITE THE I/O FOR TIME-DOMAINE DATA
load([work_dir '\VTEM_data_DF'])
load([work_dir '\VTEM_xyz_DF'])
load([work_dir '\VTEM_tc_DF'])
load([work_dir '\VTEM_Waveform'])

% Number of time channels to skip
early_tc = 0;
late_tc = 1;
radius = 13;

pT = 1e-12;
A = pi*radius^2;
% data = (data)*pT*A;

%% AEROTEM SURVEY
% load([work_dir '\ATEM_data_DF'])
% load([work_dir '\ATEM_xyz_DF'])
% load([work_dir '\ATEM_tc_DF'])
% load([work_dir '\ATEM_Waveform'])
% 
% % Number of time channels to skip
% early_tc = 0;
% late_tc = 0;
% radius = 2.5;
% 
% nT = 1e-09;
% A = pi*radius^2;
% NI = 250*8;
% data = -(data) * nT / NI;

%% Apply correction

data = data(:, (1+early_tc) : (end - late_tc) );
tc = tc((1+early_tc) : (end - late_tc));
ntc = length(tc);

std_data = std(data,1);

%% STACK DATA

%% OR DOWNSAMPLE ALONG LINE
% limits(1,1:2) = [556800 7133250];
% limits(2,1:2) = [557800 7133900];

limits(1,1:2) = [556850 7133500];
limits(2,1:2) = [557820 7134500];

[data,xyz,uncert] = sort_EM1DTM(data,xyz,limits,20,'pos');

%% Get line numbers
lineID = xy_2_lineID(xyz(:,1),xyz(:,2));
line = unique(lineID);



%% Write to H3D data format
% write_h3d_obs([work_dir '\VTEM_AVG_h3d.obs'],data,xyz,tc,radius,uncert,0.05);

%% Plot all data lines
% count = 1;
% 
% for ii = 1 : length(line)
%     set(figure(ceil(ii/4)+1), 'Position', [0 0 2000 1000])
%     subplot(4,1,count)
%     for jj = 1 : 2 : ntc
%         
%         semilogy(xyz(lineID==line(ii),1),data(lineID==line(ii),jj),:); hold on
%         
%     end
%     
%     count = count+1;
%     
%     if count==5
%         
%         count =1;
%         
%     end
%     
%     
% end


%% Create Querry (Q) matrix  for 1D to 3D mesh
% Q = make_EM1D_Q_3D(work_dir,meshfile,nullcell,xyz);
% save([work_dir '\Q'],'Q')
% load([work_dir '\Q']);

    
%% Create interpolation (P) matrix
nnodes  = 8; % Number of nearest neighbours to interpolate with
% P = make_EM1D_P_3D(work_dir,meshfile,Q,nnodes,15,100,1,2);
% save([work_dir '\P'],'P')
% load([work_dir '\P']);

fprintf('Creating Interpolation matrix\n')
fprintf('This may take a while\n')
Q = make_EM1D_Q_3D(work_dir,meshfile,nullcell,xyz);

[P,W,indx] = make_EM1D_P_3D(work_dir,meshfile,Q,interp_n,interp_r,interp_s,alpha_con(2),alpha_con(3));


%% RUN 1D INVERSION

%% Create transmiter loop (Seogi's func)
xc = 0.;
yc = 0.;
[x, y] = circfun(xc, yc, radius, 24);
txloc = [x(:), y(:)];
rxloc = [0 0];

% Input#2
% txheight = height of TX (negative above the surface)
ft = 0.3048;
txheight = xyz(:,4);
rxheight = xyz(:,4);

% Create inform background models for now


% Input#6
% wf = waveform, can be 'STEPOFF' for stepoff or a number for RAMP or
% two-column for a discretized waveform (time in second)
% ta = 5.5*1e-4;
% tb = 1.1*1e-3;
% twave = linspace(0., tb, 2^7+1) ;
% wfval = [trifun(twave, ta, tb); 0 ];

% twave = [twave 0.0026];

% waveform = [twave(:) , wfval(:) ]; 

% Input#7
% tc = vector for time channels 


%% Write EM1DTM data to file
obsfile = 'VTEM_TKC_Sub_RAW.obs';
write_EM1DTM_obs(work_dir,obsfile,-data,tc,txloc,xyz,ones(size(data)))
    
%% Setup inversion
mcell = length(nullcell);
% Pre-allocate for output model
m_out3D = ones(mcell,1) * 1e-8;

% % Create derivative matrices
% [Wx, Wy, Wz] = make_EM1DFM_grad(work_dir,meshfile);

% Pre-set average misfit
phid_avg = 99999;
beta = ones(size(Q,1),1)*1e+2;
phid = ones(size(Q,1),1)*99999;


ndata = numel(data);
% noisefun = @(a, t)(1-exp(-a*t));
% floormax = 1.5/(NI);
% floor = (min(abs(data),2))*0.025;%abs(noisefun(1000, time_chl)*floormax)*1e-2;%
% 
% sd = floor;
pred = zeros(size(data,1),size(data,2));

for ii= 1 : 7 %phid_avg(end) > target * 2;%phid_avg(end) > target

    %Leave all weights to 1 for now
    w=ones(length(nullcell),1);
%     save([work_dir '\w_iter' num2str(ii) '.con'],'-ascii','w')
    
    uncert = uncert ;
    
    % Run the inversions
    [m_con1D,d_misfit,phid,phim,beta,pred,bHSpace] = run_EM1DTM_inv(work_dir,meshfile,-data,pred,Q,phid,beta,m_con,m_ref,ii,txloc,txheight,rxloc,rxheight,uncert,'dbzdt',tc, waveform);
        
    % Interpolate between recovered column model
    m_con = P * m_con1D;

    m_ref = m_con;
    %Update model in

%     models_in{2} = m_con;
%     if ii==1
%         bHSpace = P*bHSpace;
%         
%        save([work_dir '\Bestfitting_HS.dat'],'-ascii','bHSpace');
%        
%     end
    
    % Cut topography out and save model
    m_out3D(nullcell==1) = m_con(nullcell==1);
    m_out3D(m_out3D==0) = 1e-8;
    save([work_dir '\VTEM_Inv_1D.con'],'-ascii','m_out3D')
    
    % Save misfit map
    misfit_out = P * d_misfit;
    save([work_dir '\Misfit_.con'],'-ascii','misfit_out')
    
    phid_avg(ii) = mean(phid);
    
    %% Get rig of all soundings with high misfit
    index = phid<(median(phid)*10);    
    
%     figure(2)
%     for jj = 1 : length(line)
%         
%         subplot(5,2,jj)
%         semilogy(xyz(lineID==line(jj),1),phid(lineID==line(jj)));hold on
%         
%     end
    
 
end

%% Plot obs vs pred interpolated in 2D

for ii = 1 : 5 : ntc-2
    set(figure, 'Position', [0 0 2000 1000])
    if ii == 1
        % Set coordinates for plot
        xmin = min(xyz(:,1));
        xmax = max(xyz(:,1));

        ymin = min(xyz(:,2));
        ymax = max(xyz(:,2));

        dx = 10;
        dy = 10;

        x = xmin:dx:xmax;
        y = ymin:dy:ymax;
        [Y,X] = ndgrid(y,x);

        Y = flipud(Y);
    end
    
    F_d = TriScatteredInterp(xyz(:,2),xyz(:,1),log10(data(:,ii)),'linear');
    F_p = TriScatteredInterp(xyz(:,2),xyz(:,1),log10(-pred(:,ii)),'linear');
    
    data_2D = F_d(Y,X);
    pred_2D = F_p(Y,X);
    
    subplot(1,3,1)
    imagesc(y,x,data_2D);
    title(['\bf LOG Observed Time: ' num2str(tc(ii))]);
%     caxis([min(data(:,ii)) max(data(:,ii))]);
    colorbar
    
    subplot(1,3,2)
    imagesc(y,x,pred_2D);
    title(['\bf LOG Predicted : ' num2str(tc(ii))]);
%     caxis([min(data(:,ii)) max(data(:,ii))]);
    colorbar
    
    subplot(1,3,3)
    imagesc(y,x,data_2D-pred_2D);
    title(['\bf LOG Residual Time: ' num2str(tc(ii))]);
%     caxis([min(data(:,ii)) max(data(:,ii))]);
    colorbar
    
end   

%% Plot obs vs pred for station inversion 
count = 1;
% 
% for ii = 1 : length(line)
%     set(figure(ceil(ii/4)+1), 'Position', [0 0 2000 1000])
%     subplot(4,1,count)
%     for jj = 1 : 2 : ntc
%         
%         semilogy(xyz(lineID==line(ii),1),data(lineID==line(ii),jj),':'); hold on
%         semilogy(xyz(lineID==line(ii),1),-pred(lineID==line(ii),jj),'r'); hold on
%         
%     end
%     ylim([1e-14 1e-8])
%     count = count+1;
%     title(['\bfNorthing: ' num2str( mean( xyz(lineID==line(ii),2) ) )]);
%     if count==5
%         
%         count =1;
%         
%     end
%     
%     
% end

% set(figure, 'Position', [25 100 1800 900])
% for ii = 1 : size(data,1)
%     
%     
%     subplot(size(stn,1)/2,2,ii)
%     loglog(tc,abs(data(ii,:)),'b'); hold on
%     loglog(tc,abs(pred(ii,:)),'r'); hold on
%     legend('Obs','Pred')
%     title(['\bfEasting:' num2str(stn(ii,1)) 'm Northing:' num2str(stn(ii,2)) 'm']) 
%     xlim([min(tc) max(tc)])
%     xlabel('\bfTime (s)')
%     ylabel('\bfdB/dt (nT/s)')
%     grid on
%     
% end
% axis([min(tc) max(tc) abs(max(data(:,1))) abs(min(data(:,end)))])
