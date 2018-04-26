% Driver for DHEM

clear all
close all


root_dir = pwd;
addpath functions
addpath 'C:\Users\DominiqueFournier\Documents\GIT\ScriptsBank\FUNC_LIB'

% Project path
basePath = 'C:\Users\DominiqueFournier\Google Drive\Research\Hudbay_Lalor';
fwr_dir = [basePath '\Modeling\Forward\SQUID'];
timefile = 'times_out.txt';
FWRfile = 'recv_h3dtd.txt';
predfile = 'Modeling\Inversion\Inv20_avgdata\dpred_10.txt';

% work_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\PEM_test';
work_dir = [basePath '\Processing'];

result_dir = [basePath  '\Modeling\Inversion\Inv20_avgdata'];

% data_dir = 'C:\LC\Private\dominiquef\Projects\3796_AGIC_Research\PEM_test\PEM_raw';
data_dir = [basePath  '\Data\Crone_data\Phase2'];

meshfile = [basePath  '\Inversion\Inv3_fine_noredata\DHEM_mesh_fine.msh'];

% Extract information from PEM files
DHEM_raw = Read_PEM(data_dir);
save([work_dir '\' 'Lalor_DHEM_RAW'],'DHEM_raw');
load ([work_dir '\' 'Lalor_DHEM_RAW'])

% De-survey each hole and interpolate between stations for dip and azimuth
Atlas = DAD2Atlas(work_dir,DHEM_raw);
save([work_dir '\' 'Lalor_DHEM_UTM'],'Atlas');
load ([work_dir '\' 'Lalor_DHEM_UTM'])


%% INPUT PARAMETERS
nstn = size(Atlas.dBx,1);
ndh = size(DHEM_raw.name,1);

% t0 = 1.5/1000;
% first_time = 2;
% last_time = 25;
% subtime = [11 13 15 17 19 20 21 22 23]';

t0 = 0.1515;
first_time = 3;
last_time = 37;
subtime = [1:34]';

%% WRITE TRANSMITTER TO FILE (TEST)
% fid = fopen([work_dir '\' 'Tx_loop.dat'],'w');
% 
% for jj = 1:size(DHEM_UTM.tx{1},1)
%     fprintf(fid,'VRTX %i %12.3f %12.3f %12.3f\n',jj, DHEM_UTM.tx{1}(jj,1),...
%         DHEM_UTM.tx{1}(jj,2),DHEM_UTM.tx{1}(jj,3));
% end
% 
% for jj = 1:size(DHEM_UTM.tx{1},1)-1
%     fprintf(fid,'SEG %i %i\n',jj, jj+1);
% end
% 
% fclose(fid);

%% WRITE TIME CHANNEL FILE

fid = fopen([work_dir '\times_out.txt'],'w');

t_channel = zeros(last_time - first_time ,1);

ntc = length(t_channel);

for ii = 1 : ntc
    
    t_channel(ii) = (DHEM_raw.time{1}(first_time+ii-1)+DHEM_raw.time{1}(first_time+ii)) / 2 + t0;
    
    fprintf(fid,'%8.6f\n',t_channel(ii));
    
end

fclose(fid);


load([fwr_dir '\wave.dat']);

%% BUILD WAVEFORM
% tt = 0.000:0.00005:0.015;
% rampup = 1+(log(100*(tt)/max(tt)));
% rampdw = @(t) exp(-(t)/0.000125);
% 
% figure; plot([1e-3 0.1509 tt+0.1509],[1 1 rampdw(tt)]);ylim([0 1.1]); xlim([0.14 0.3]);hold on
% plot(t_channel,rampdw(t_channel-0.015),'g*')
% hold on;
% plot(wave(:,1),wave(:,2),'ro');xlim([0.15075 0.153])

%% Build tx file for forward model
% 
% fid = fopen([work_dir '\trx_loc.txt'],'w');
% fprintf(fid,'N_TRX    1\n');
% fprintf(fid,'\nTRX_ORIG\n');
% 
% ntx = size(DHEM_UTM.tx{1},1);
% fprintf(fid,'%i\n',ntx);
% 
% for ii = 1 : ntx
%     
%     fprintf(fid,'%12.3f %12.3f %12.3f\n',DHEM_raw.tx{1}(ii,1),DHEM_raw.tx{1}(ii,2),DHEM_raw.tx{1}(ii,3));
%     
% end
% 
% fprintf(fid,'INCLUDE recv_loc.txt\n');
% fclose(fid);


%% Build mesh for interp
% [xn,yn,zn] = read_UBC_mesh(meshfile)
% 
% nx = mesh(1,1);
% ny = mesh(1,2);
% nz = mesh(1,3);
% 
% ncell = nx*ny*nz;
% 
% dx = mesh(3,1:nx);
% dy = mesh(4,1:ny);
% dz = mesh(5,1:nz);
% 
% % Create cell center coordinate vector
% x = mesh(2,1) + cumsum(dx) - dx/2;
% y = mesh(2,2) + cumsum(dy) - dy/2;
% z = mesh(2,3) - cumsum(dz) + dz/2;
% 
% [Z,X,Y] = ndgrid(z,x,y);

xmin = min(Atlas.DXYZ(:,2));
xmax = max(Atlas.DXYZ(:,2));

ymin = min(Atlas.DXYZ(:,3));
ymax = max(Atlas.DXYZ(:,3));

zmin = min(Atlas.DXYZ(:,4));
zmax = max(Atlas.DXYZ(:,4));



%% Downsample the data for inversion - By distance along hole
% or by grid cell size
index = ones(size(Atlas.DXYZ,1),1);

radius = 5;

dobs = zeros(1,22);
recv = zeros(1);

dHid = [];
count = 1;
count_rx = 0;

global_std_dBx = std(Atlas.dBx);
global_std_dBy = std(Atlas.dBy);
global_std_dBz = std(Atlas.dBz);

for ii = 1 : size(Atlas.DXYZ,1)
           
    if index(ii)==1
        
        count_rx = count_rx + 1;
%         dobs(count,:) = NaN(1,22);
        
%         temp = zeros(length(rc_specs),1); temp(ii:ii+nfreq-1)=1;
%             temp = stn_num(:,1)==stn(ii,1);
        r = ( ( Atlas.DXYZ(ii,2) - Atlas.DXYZ(:,2)).^2 +...
            (Atlas.DXYZ(ii,3) - Atlas.DXYZ(:,3)).^2 +...
            (Atlas.DXYZ(ii,4) - Atlas.DXYZ(:,4)).^2 ) .^0.5;

        % Only keep the closest to cell center
        index(r <= radius) = 0;
        index(r==(min(r))) = 1;


        for jj = 1 : ntc
        dobs(count,1:4) = [Atlas.DXYZ(ii,2) Atlas.DXYZ(ii,3) Atlas.DXYZ(ii,4) t_channel(jj)];
        
        dobs(count,17:22) = [Atlas.dBx(ii,jj) global_std_dBx(jj) ...
                            Atlas.dBy(ii,jj) global_std_dBy(jj) ...
                            Atlas.dBz(ii,jj) global_std_dBz(jj) ];
                        
        temp = unique(Atlas.Htag(ii));
        dHid(count) = temp(1);               
        recv(count) = count_rx;
        
        count = count + 1;
                     
        end
        
    end
  
end

dobs_flt = dobs ;
dHid_flt = dHid ;
save([work_dir '\Subset_index'],'index');
save([work_dir '\Obs_filter_25m_DH'],'dobs_flt');
save([work_dir '\holeID_flt_25m'],'dHid_flt');
save([work_dir '\Rx_25m'],'recv');

load([work_dir '\Subset_index']);
load([work_dir '\holeID_flt_25m']);
load([work_dir '\Obs_filter_25m_DH']);
load([work_dir '\Rx_25m']);

%% DOWN-SAMPLE DATA
% Select data within cells and compute weighted average
% 
% % Create global standard deviation in case if not enough obs within a cell
% global_std_dBx = std(Atlas.dBx);
% global_std_dBy = std(Atlas.dBy);
% global_std_dBz = std(Atlas.dBz);
% 
% % Define interpolation grid
% dx = 50;
% dy = 50;
% dz = 50;
% 
% x = floor(xmin-dx/2) : dx : ceil(xmax+dx/2); nx = length(x);
% y = floor(ymin-dy/2) : dy : ceil(ymax+dy/2); ny = length(y);
% z = ceil(zmax+dz/2) : -dz : floor(zmin-dz/2); nz = length(z);
% 
% dx = ones(nx,1)*dx;
% dy = ones(ny,1)*dy;
% dz = ones(nz,1)*dz;
% 
% [Z,X,Y] = ndgrid(z,x,y);
% 
% % Pre-allocate for the obs, but will grow as more average data is added
% dobs = zeros(1,22);
% dHid = [];
% count = 1;
% nrecv = 0;
% for kk = 1:ny
%     
%     for jj = 1:nx
%         
%         for ii = 1:nz
%             
%             
%             
%             indexx = ((X(ii,jj,kk) - dx(jj)/2) <= Atlas.DXYZ(:,2)) .* ((X(ii,jj,kk) + dx(jj)/2) >= Atlas.DXYZ(:,2));
%             indexy = ((Y(ii,jj,kk) - dy(kk)/2) <= Atlas.DXYZ(:,3)) .* ((Y(ii,jj,kk) + dy(kk)/2) >= Atlas.DXYZ(:,3));
%             indexz = ((Z(ii,jj,kk) - dz(ii)/2) <= Atlas.DXYZ(:,4)) .* ((Z(ii,jj,kk) + dz(ii)/2) >= Atlas.DXYZ(:,4));
%             
%             index = indexx .* indexy .* indexz == 1;
%             
%             if sum(index) > 1
%                 
%                 nrecv = nrecv + 1 ;
%                 % Extract data around location
%                 sub_dBx = Atlas.dBx(index,:);
%                 sub_dBy = Atlas.dBy(index,:);
%                 sub_dBz = Atlas.dBz(index,:);
% 
%                 % Pre-allocate matrix
%                 avg_dBx = zeros(1,ntc);
%                 avg_dBy = zeros(1,ntc);
%                 avg_dBz = zeros(1,ntc);
% 
%                 % Computed weighted average
%                 for tt = 1 : ntc
% 
%                     dobs(count,:) = NaN(1,22);
%                  % Only keep positive data for each time channel
%                 %          stck(:,jj) = abs(stck(:,jj));
%                 %          select = stck(:,jj)>0;
% 
%                     std_dBx = std(sub_dBx(:,tt));
%                     std_dBy = std(sub_dBy(:,tt));
%                     std_dBz = std(sub_dBz(:,tt));
% 
%                     mux = mean(sub_dBx(:,tt));
%                     muy = mean(sub_dBy(:,tt));
%                     muz = mean(sub_dBz(:,tt));
% 
%                     wx = 1./abs(sub_dBx(:,tt) - mux );
%                     wy = 1./abs(sub_dBy(:,tt) - muy );
%                     wz = 1./abs(sub_dBz(:,tt) - muz );
% 
%                     avg_dBx = sum(wx.*sub_dBx(:,tt)) / sum(wx);
%                     avg_dBy = sum(wy.*sub_dBy(:,tt)) / sum(wy);
%                     avg_dBz = sum(wz.*sub_dBz(:,tt)) / sum(wz);
% 
%                 %              semilogx(tc(jj),mu,'ro'); hold on
%                 %              semilogx(tc(jj),stck(:,jj),'*','MarkerSize',2); hold on
%                     dobs(count,1:4) = [X(ii,jj,kk) Y(ii,jj,kk) Z(ii,jj,kk) t_channel(tt)];
%                     
%                     if sum(index) > 5
%                         dobs(count,17:22) = [avg_dBx std_dBx ...
%                                              avg_dBy std_dBy ...
%                                              avg_dBz std_dBz ];
% 
%                           temp = unique(Atlas.Htag(index));
%                         dHid(count) = temp(1);
% 
%                     else
%                         
%                         dobs(count,17:22) = [avg_dBx global_std_dBx(tt) ...
%                                              avg_dBy global_std_dBy(tt) ...
%                                              avg_dBz global_std_dBz(tt) ];
%                            temp = unique(Atlas.Htag(index));
%                         dHid(count) = temp(1);              
%                     end
% 
%                      count = count + 1;
%                 end
% 
%             end
%             
%         end
%         
%     end
%     
% end
% 
% 
% save([work_dir '\Obs_avg_40m'],'dobs');
% save([work_dir '\holeID'],'dHid');
% load([work_dir '\holeID']);
% load([work_dir '\Obs_avg_40m']);

%% OR TAKE ALL DATA

% ntc = length(t_channel);
% nrecv = size(Atlas.dBx,1);
% 
% dobs = ones(ntc*nrecv,22)*-99999;
% 
% count = 0;
% for ii = 1 : nrecv
%     
%     for jj = 1 : ntc
%         count = count+1;
%         
%         dobs(count,:) = [Atlas.DXYZ(ii,2),Atlas.DXYZ(ii,3),Atlas.DXYZ(ii,4),t_channel(jj)...
%             (ones(1,12)*-99),Atlas.dBx(ii,jj),abs(Atlas.dBx(ii,jj))*0.25...
%             Atlas.dBy(ii,jj),abs(Atlas.dBy(ii,jj))*0.25...
%             Atlas.dBz(ii,jj),abs(Atlas.dBz(ii,jj))*0.25];
%         
%     end
%     
% end
% 
% dobs_flt = dobs;
% % Create receiver ID array
% nrecv = size(dobs,1) / ntc;
% recv = kron([1:nrecv]',ones(ntc,1));
% index = ones(nrecv,1);

%% WRITE OBS FILE
ndv = 99999;
pc_err = 0.025;

ntc = length(t_channel);

% Select time channels 
time_out = t_channel(subtime) ;

% Components [1=x,2=y,3=z]
component = [1 2 3];

% Transmitter loop
tx = DHEM_raw.tx{1};



% Write rec_loc.dat
index = index==1;
recv_xyz = [Atlas.DXYZ(index,2) Atlas.DXYZ(index,3) Atlas.DXYZ(index,4)];
save([work_dir '\recv_loc.dat'],'-ascii','recv_xyz');
save([work_dir '\times_out.txt'],'-ascii','time_out');

write_DHEM_obs(work_dir,'DHEM_obs.dat',dobs_flt,tx,recv',subtime,component,pc_err)



%% ! ! RUN FORWARD MODEL OR INVERSION NOW ! !

%% Load forward model from file and PLOT
tc = load([work_dir '\' timefile]);
ntc_out = length(tc);

%% LOAD obs and interpolate in 3D
% 
% dx = 20;
% dy = 20;
% dz = 20;
% 
% x = floor(xmin-dx/2) : dx : ceil(xmax+dx/2); nx = length(x);
% y = floor(ymin-dy/2) : dy : ceil(ymax+dy/2); ny = length(y);
% z = ceil(zmax+dz/2) : -dz : floor(zmin-dz/2); nz = length(z);
% 
% dx = ones(nx,1)*dx;
% dy = ones(ny,1)*dy;
% dz = ones(nz,1)*dz;
% 
% write_UBC_mesh(work_dir,xmin,ymin,zmax,...
%     dx,dy,dz);
% 
% [Z,X,Y] = ndgrid(z,x,y);
% 
% mcell = nx * ny * nz;
% 
% obsfile = 'DHEM_obs_avg50m.dat';
% 
% % Extract complete obs and pred and write to Gocad file if argin == write
% index = 1:nstn;
% index = find(index>=0);
% % [dhid,dobs,dpred,dstd,dres,dloc] = H3D_Pred_vs_obs(Atlas,fwr_dir,obsfile,predfile,index,'none');
% 
% % Read obsfile
% % [dobs] = read_H3D_obs([result_dir '\' obsfile]);
% 
% % Pre-allocation
% Interp_dBxdt = zeros(mcell,ntc_out);
% Interp_dBydt = zeros(mcell,ntc_out);
% Interp_dBzdt = zeros(mcell,ntc_out);
% 
% % Interpolate all time channels onto 3D mesh
% for ii = 1 : 10:11 %ntc
%     
%     
%     xobs = dobs(ii:ntc:end,1);
%     yobs = dobs(ii:ntc:end,2);
%     zobs = dobs(ii:ntc:end,3);
%     obsx = dobs(ii:ntc:end,17);
%     obsy = dobs(ii:ntc:end,19);
%     obsz = dobs(ii:ntc:end,21);
%     
%     fprintf('Interpolating Obs dBx/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsx,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %     temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     
%     Interp_dBxdt(:,ii) = temp(:);
%     save([result_dir '\ObsX_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
%     fprintf('Interpolating Obs dBy/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsy,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %     temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     
%     Interp_dBydt(:,ii) = temp(:);
%     save([result_dir '\ObsY_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
%     fprintf('Interpolating Obs dBz/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsz,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %         temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     
%     Interp_dBzdt(:,ii) = temp(:);
%     save([result_dir '\ObsZ_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
% end

%% Load predicted and interpolate on 3D mesh

% predfile = 'dpred_10.txt';
% ntc_out = length(time_out);
% % Extract complete obs and pred and write to Gocad file if argin == write
% % Read pred
% dpre = load([result_dir '\' predfile]);
% 
% % Pre-allocation
% Interp_dBxdt = zeros(mcell,ntc_out);
% Interp_dBydt = zeros(mcell,ntc_out);
% Interp_dBzdt = zeros(mcell,ntc_out);
% 
% % Interpolate all time channels onto 3D mesh
% for ii = [1 ntc_out] %ntc
%     
%     fprintf('Interpolating PRED at time: %f\n',tc(ii));
%     xobs = dpre(ii:ntc_out:end,1);
%     yobs = dpre(ii:ntc_out:end,2);
%     zobs = dpre(ii:ntc_out:end,3);
%     obsx = dpre(ii:ntc_out:end,11);
%     obsy = dpre(ii:ntc_out:end,12);
%     obsz = dpre(ii:ntc_out:end,13);
%     
%     fprintf('Interpolating PRED dBx/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsx,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %     temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     Interp_dBxdt(:,ii) = temp(:);
%     save([result_dir '\Inv20_predX_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
%     fprintf('Interpolating PRED dBy/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsy,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %     temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     Interp_dBydt(:,ii) = temp(:);
%     save([result_dir '\Inv20_predY_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
%     fprintf('Interpolating PRED dBz/dt at time: %f\n',tc(ii));
%     temp = griddata(xobs,yobs,zobs,obsz,X(:),Y(:),Z(:),'natural');
% %     temp = log10( abs( temp) );
% %     temp(temp > 0) = temp(temp > 0)/max(max(max(temp(temp > 0))));
% %     temp(temp < 0) = -temp(temp < 0)/min(min(min(temp(temp < 0))));
%     temp(isnan(temp)) = -99999;
%     Interp_dBzdt(:,ii) = temp(:);
%     save([result_dir '\Inv20_predZ_tc' num2str(ii) '.dat'],'-ascii','temp');
%     
% end



%% RUN FORWARD MODEL
%% LOAD DATA AND PLOT
% dobs = read_H3D_obs([work_dir '\DHEM_obs_v2.dat']);
dpre = load([fwr_dir '\' FWRfile]);
% tc = load([work_dir '\times_out.txt']);
close all
plot_hole = [1]%randi(length(unique(dHid)),1,5);

for ii = 1 : length(plot_hole)
        
        
    for jj = 1 : length(time_out);
    set(figure(jj), 'Position', [50 50 500 800]);
    t_obs = dobs(:,4) == time_out(jj);
    
    
    F_obsx = TriScatteredInterp(dobs(t_obs,1),dobs(t_obs,2),dobs(t_obs,3),dobs(t_obs,17));
    F_obsy = TriScatteredInterp(dobs(t_obs,1),dobs(t_obs,2),dobs(t_obs,3),dobs(t_obs,19));
    F_obsz = TriScatteredInterp(dobs(t_obs,1),dobs(t_obs,2),dobs(t_obs,3),dobs(t_obs,21));
    
        
    t_pre = round(dpre(:,4)*1e+5) == round(time_out(jj)*1e+5);
    
    F_prex = TriScatteredInterp(dpre(t_pre,1),dpre(t_pre,2),dpre(t_pre,3),dpre(t_pre,11));
    F_prey = TriScatteredInterp(dpre(t_pre,1),dpre(t_pre,2),dpre(t_pre,3),dpre(t_pre,12));
    F_prez = TriScatteredInterp(dpre(t_pre,1),dpre(t_pre,2),dpre(t_pre,3),dpre(t_pre,13));
    count = 1;    
    
            
        holeid = plot_hole(ii);
%         grab1 = dHid'==holeid;
        grab2 = dHid_flt'==holeid;
    
        
%         subplot(1,3,ii)

%         dobs_1 = dobs( dobs(:,4) == t_channel(jj) & grab1 , : );
        obs_DH = dobs_flt( dobs_flt(:,4) == time_out(jj) & grab2 , : );
        
        
        obs_dBx = F_obsx(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        obs_dBy = F_obsy(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        obs_dBz = F_obsz(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        
        pre_dBx = F_prex(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        pre_dBy = F_prey(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        pre_dBz = F_prez(obs_DH(:,1),obs_DH(:,2),obs_DH(:,3));
        
%         dpre_t = dpre( dpre(:,4) == dpre(jj,4), :);

        plot( obs_dBx , obs_DH(:,3) - 257 , 'r:','LineWidth',1) ; hold on     
        plot( obs_dBy , obs_DH(:,3) - 257 , 'g:','LineWidth',1) ; hold on        
        plot( obs_dBz , obs_DH(:,3) - 257 ,  'b:','LineWidth',1) ; hold on
        
        plot( -pre_dBx , obs_DH(:,3) - 257 , 'r','LineWidth',2) ; hold on     
        plot( -pre_dBy , obs_DH(:,3) - 257 , 'g','LineWidth',2) ; hold on        
        plot( -pre_dBz , obs_DH(:,3) - 257 , 'b','LineWidth',2) ; hold on
        grid on
        ylim([-1000,25]);
%         plot( obs_DH(:,19) , obs_DH(:,3) , 'g:','LineWidth',2) ; hold on
%         plot( obs_DH(:,21) , obs_DH(:,3) ,  'b:','LineWidth',2) ; hold on
%         plot( obs_DH(:,17) , obs_DH(:,3) ,  'r:','LineWidth',2) ; hold on

%         plot( (dpre_t(grab,11)), dpre_t(grab,3)  , ':r') ; hold on
%         plot( (dpre_t(grab,12)), dpre_t(grab,3)  , ':g') ; hold on
%         plot( (dpre_t(grab,13)), dpre_t(grab,3)  , ':b') ; hold on

%         subplot(1,3,ii)
        title(['\bf Hole:' DHEM_raw.name{plot_hole(ii)}])
        xlabel('\bfAmplitude (T/s)')
        ylabel('\bfDepth (m)')
        legend('Obs dBx/dt','Obs dBy/dt','Obs dBz/dt','Pred dBx/dt','Pred dBy/dt','Pred dBz/dt','Location','SouthEast')
        count = count+1;
    end

end

%% Plot data in slice and forward fields
% Create logical vector for remaining data stations
% in = keep1==1;
% % slices = [6 10 14 18];
% % 
% % y = y(7:end-7);
% % x = x(10:end-15);
% % 
% % [Z,X,Y] = ndgrid(z,x,y);
% % 
% % [F] = DHEM_plot(FWR_times,FWR,Atlas,X,Y,Z,in,slices,'natural');
% % save([result_dir '\' 'SQUID'],'F')
% 
% load([result_dir '\' 'SQUID']) 
% [h, w, p] = size(F(1).cdata);
% hf = figure;
% 
% set(hf,'Position', [25 100 1800 900]);
% axis off
% 
% movie(hf,F,20,10,[0 0 0 0])


%% Plot slices of interpolated data
% x = floor(xmin-dx/2) : 10 : ceil(xmax+dx/2); nx = length(x);
% y = floor(ymin-dy/2) : 10 : ceil(ymax+dy/2); ny = length(y);
% [X,Y] = ndgrid(x,y);
% count = 1;
% ntc_sub = length(time_out);
% 
% % Load pred
% % dpre = load('C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Hudbay_Lalor\Processing\FWR_SQUID\recv_h3dtd.txt');
% count = 0;
% for ii = 11
%     
%     subdata = dobs(ii:ntc_sub:end,:);
%     subpre = dpre(ii:ntc_sub:end,:);
% 
% %     
%     for jj = 15 : 5 : 25
%         
%         count = count + 1;
%         index = subdata(:,3) == Z(jj,1,1);
%         
%         if sum(index)~= 0
%             
%         set(figure(count), 'Position', [00 00 2000 1800]);
%         subplot(2,3,1)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subdata(index,17),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         
%         imagesc(x,y,d_slice);
%         caxis([-1 1]) 
%         colorbar
%         title(['\bfOBS dBdx T:' num2str(time_out(ii)) 's'])
%        
%         subplot(2,3,2)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subdata(index,19),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         imagesc(x,y,d_slice);
%       caxis([-1 1]) 
%         colorbar
%         title(['\bfOBS dBdy T:' num2str(time_out(ii)) 's'])
%         
%         subplot(2,3,3)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subdata(index,21),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         imagesc(x,y,d_slice);
%         caxis([-1 1]) 
%         colorbar
%         title(['\bfOBS dBdz T:' num2str(time_out(ii)) 's'])
%         
%         subplot(2,3,4)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subpre(index,11),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         imagesc(x,y,d_slice);
%         caxis([-1 1]) 
%         colorbar
%         title(['\bfPRE dBdx T:' num2str(time_out(ii)) 's'])
%         
%         subplot(2,3,5)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subpre(index,12),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         imagesc(x,y,d_slice);
%         caxis([-1 1]) 
%         colorbar
%         title(['\bfPRE dBdy T:' num2str(time_out(ii)) 's'])
%         
%         subplot(2,3,6)
%         d_slice = griddata(subdata(index,1),subdata(index,2),subpre(index,13),X,Y,'natural');
%         d_slice(d_slice > 0) = d_slice(d_slice > 0)/max(d_slice(d_slice > 0));
%         d_slice(d_slice < 0) = -d_slice(d_slice < 0)/min(d_slice(d_slice < 0));
%         imagesc(x,y,d_slice);
%         caxis([-1 1]) 
%         colorbar
%         title(['\bfPRE dBdz T:' num2str(time_out(ii)) 's'])
%         end
%         
%     end
%     
%     count = count + 1;
%     
% end
        
        


