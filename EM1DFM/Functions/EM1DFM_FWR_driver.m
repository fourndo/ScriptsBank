% FOR DEV ONLY
% Temporary function to format data to GIFTool format
% MUST BE REFORMATED FOR EVERY FILE

clear all
close all

addpath     'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB';

work_dir    = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Tli_Kwi_Cho\Modelling\Inversion\EM\DIGHEM\1D\FWR_VTEM';
meshfile    = 'VTEM1D_mesh_small.msh';
confile    = 'VTEM_Inv_1D.con';
obsfile     = 'DIGHEM_TKC_ALL.obs';
topofile    = 'CDED_076c05_NAD27.topo';

%% Reformat data in array structure
% Last entry specifies the minimum distance between points
freqin = [56000 7200 900];% 5001 901];
limits(1,1:2) = [556800 7133250];
limits(2,1:2) = [557800 7133900];

[data,xyz] = rawdata_2_EM1DFM([work_dir '\DIGHEM_data'],freqin,25,limits);


% Remove all zero observation
% ndv = data{7}(:,1)~=0 & data{7}(:,2)~=0;
% for ii = 1 : size(data,2)
%     
%     data{ii} = data{ii}(ndv,:);
%     
% end

%% Write EM1DFME data to file

writeem1dfmobs(work_dir,obsfile,data,'')

% %% Write lines of data
% % Assign line number to data
lineID = xy_2_lineID(xyz(:,1),xyz(:,2));
line = unique(lineID);

% Take a subset

% for ii = 1 : length(line);
%     
%     for jj = 1 : size(data,2)
%         
%         subdata{jj} = data{jj}(lineID==line(ii),:);
%         
%     end
%     
%     writeem1dfmobs(work_dir,['DIGHEM_line' num2str(line(ii)) '.obs'],subdata,'')
% 
% end

%% Make nullcell from topo and mesh
% [xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);
% [Zn,Xn,Yn] = ndgrid(zn,xn,yn);
% 
% % Load topography
% topo = read_UBC_topo([work_dir '\' topofile]);
% 
% % Create discretize topogrphy
% [nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% save([work_dir '\nullcell.dat'],'-ascii','nullcell');
load([work_dir '\nullcell.dat']);
% m_sus_ref = load([work_dir '\' sus_file]);

%% Create Querry (Q) matrix  for 1D to 3D mesh
Q = make_EM1D_Q_3D(work_dir,meshfile,nullcell,xyz);
save([work_dir '\Q'],'Q')
load([work_dir '\Q']);


%% RUN 1D INVERSION

mtype = '1';
% Create inform background models for now
m_con = load([work_dir '\' confile]); 
% m_con(m_con==1) = 1e-4;
% m_con(m_con==0) = 1e-4;

m_sus = nullcell; 
m_sus(m_sus==1) = 0;

% m_sus_ref = m_sus;

mcell = length(nullcell);
 
% Run the inversions
[pred] = run_EM1DFM_fwr(work_dir,meshfile,obsfile,m_con,m_sus,Q);
        
lineID = kron(lineID,ones(3,1));  
%% Plot obs vs predicted



for ii = 1 : length(line)

    index = lineID == line(ii);

    set(figure(ii), 'Position', [25 50 1800 900]) 
    subplot(211)
    plot(data{9}(index & data{3}(:,1)==56000, 1), data{7}(index & data{3}(:,1)==56000, 1) , 'r:'); hold on
    plot(data{9}(index & data{3}(:,1)==7200, 1), data{7}(index & data{3}(:,1)==7200, 1) , 'b:'); hold on
    plot(data{9}(index & data{3}(:,1)==900, 1), data{7}(index & data{3}(:,1)==900, 1) , 'g:'); hold on
    
    plot(pred(index & pred(:,4)==56000, 1), pred(index & pred(:,4)==56000, 5) , 'r'); hold on
    plot(pred(index & pred(:,4)==7200, 1), pred(index & pred(:,4)==7200, 5) , 'b'); hold on
    plot(pred(index & pred(:,4)==900, 1), pred(index & pred(:,4)==900, 5) , 'g'); hold on
    
    subplot(212)
    plot(data{9}(index & data{3}(:,1)==56000 ,1), data{7}(index & data{3}(:,1)==56000, 2) , 'r:'); hold on
    plot(data{9}(index & data{3}(:,1)==7200 ,1), data{7}(index & data{3}(:,1)==7200, 2) , 'b:'); hold on
    plot(data{9}(index & data{3}(:,1)==900 ,1), data{7}(index & data{3}(:,1)==900, 2) , 'g:'); hold on

    plot(pred(index & pred(:,4)==56000, 1), pred(index & pred(:,4)==56000, 6) , 'r'); hold on
    plot(pred(index & pred(:,4)==7200, 1), pred(index & pred(:,4)==7200, 6) , 'b'); hold on
    plot(pred(index & pred(:,4)==900, 1), pred(index & pred(:,4)==900, 6) , 'g'); hold on
    
end