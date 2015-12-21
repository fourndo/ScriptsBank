% Advanced Magnetic Inversion (AMI3D)
% Inversion algorithm for magnetostatic problem
%
% Starting with TMI data, computes amplitude data using Equivalent Source
% code.
%
% Joint code using Magnetic Vector Inversion (MVI) and Magnetic Amplitude
% Inversion (MAI) codes.
% 
% INPUT PARAMETER
%
% 
% DEVELOPMENT CODE
% Other sub-functions are required to run the code
% Written by: Dominique Fournier 
% Last update: May 31th, 2015

% clear all
% close all

addpath ..\FUNC_LIB\;
addpath AMI;

% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\AMI';
inpfile   = 'MAGINV_TMVI.inp'; 
topofile = []; % TO BE ADDED TO THE INPUTS

[meshfile,obsfile,mstart,mref,m_esus,chi_target,alphas,beta,bounds,lp_vec,lp_tresh,weightfile,FLAG1,FLAG2] = AMI3D_read_inp([work_dir '\' inpfile]);

% mtrue = load([work_dir '\..\Effec_sus.sus']);
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

mcell = (length(xn)-1) * (length(yn)-1) * (length(zn)-1);

pct_cutoff = 70;


%% Generate s and t vectors for the different lp zones
if ischar(lp_vec)==1
    
    lpmat = load([work_dir '\' lp_vec]);
    [s,LP] = find_zones(lpmat);

    % Smooth out the regions with 8-point averager
    % Power determine the transition length
    A = get_AVG_8pt(dx,dy,dz);
    A = A*(A*A);
    A = spdiags(1./sum(A,2),0,mcell,mcell) *A;
    
    t = A*s;
    
    
else
    
    LP = lp_vec;
    t= ones(mcell,1);

end

%% Load weights model
if isempty(weightfile) == 1
    
    w = ones(4*mcell,1);
    
else
    
    w = load([work_dir '\' weightfile]);
    
end

% Create or load reference model
if ischar(m_esus)==1
    
    m_esus = load([work_dir '\' m_esus]);
    
else
    
    m_esus = ones(mcell,1)*m_esus;
   
    
end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir '\' mstart]);
    
else
    
    mstart = ones(mcell,1)*mstart;
   
    
end

% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir '\' mref]);
    
else
    
    mref = zeros(mcell,1);
   
    
end
%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndata = length(d);
Wd   = spdiags(1./wd,0,ndata,ndata);

%% RUN EQUIVALENT SOURCE CODE
% [EMS_sus,mvec] =EMS3D(work_dir,meshfile,obsfile,topofile,mstart,mref,chi_target,alphas,beta,1,FLAG1);

EMS_sus = load([work_dir '\EMS_3D.sus']);
% mvec = load([work_dir '\EMS_mvec.fld']);


%% COOPERATIVE INVERTION
betaMAI = [];
betaMVI = [];
m = [zeros(mcell,1) zeros(mcell,1) ones(mcell,1)];

% Pre-initiate phid
phidMAI = ndata*2;
phidMVI = ndata*2;
MOF_MAI = speye(mcell);
MOF_MVI = speye(3*mcell);
esus_MAI = m_esus;
switchMAI = 3;
switchMVI = 0;
ii = 0;
beta_tol = 0.25;
while switchMAI~=3 ||  switchMVI~=3
    
    
%% RUN AMPLITUDE INVERSION (MAI)

ii = ii + 1;
fprintf('BEGIN AMI ITERATION # %i\n',ii);

    if switchMAI == 1

        max_iterMAI = 30;
        fprintf('BEGIN IRLS FOR MAI\n') 
    else

        max_iterMAI = 30;
        
    end

if switchMAI~=3
    fprintf('\n### MAI STEP\n')
%     [esus_MAI,pred_ampB,betaMAI,phidMAI,MOF_MAI,switchMAI] = MAI3D(work_dir,meshfile,'EMS_amp.pre',m_esus,mref,m,w,chi_target,phidMAI,MOF_MAI,alphas,betaMAI,[0 1],[2 2 2 2 1],t,switchMAI,FLAG1,FLAG2,max_iterMAI,beta_tol );

    
%     save([work_dir '\MAGINV_esus_iter_.sus'],'-ascii','esus_MAI')
%     write_MAG3D_TMI([work_dir '\MAG3D_esus_iter_.pre'],H,I,Dazm,obsx,obsy,obsz,pred_ampB,wd);
end
esus_MAI = load([work_dir '\MAGINV_esus_iter_.sus']);
%% PASS INFORMATION FROM MAI TO MVI
% MAI assumes all in the direction of magnetization.
[P,S,T] = azmdip_2_pst(Dazm,I,1);
Esus = spdiags(m_esus,0,mcell,mcell);

Mp = Esus * (P'*m')' ;
Ms = Esus * (S'*m')' ;
Mt = Esus * (T'*m')' ;

mstart = [Mp;Ms;Mt];
% mref = [Mp;Ms;Mt];

%% RUN MAGNETIC VECTOR INVERSION (MVI)
 
    if switchMVI == 1 && switchMAI == 3

        max_iterMVI = 30;
        fprintf('BEGIN IRLS FOR MVI\n')

    else

        max_iterMVI = 30;
        fprintf('\nMVI STEP\n')
    end

[M,pred_TMI,betaMVI,phidMVI,MOF_MVI,switchMVI] = MVI3D(work_dir,meshfile,obsfile,mstart,mref,esus_MAI,chi_target,phidMVI,MOF_MVI,alphas,betaMVI,bounds,LP,t,lp_tresh{1},lp_tresh{2} ,switchMVI,FLAG1,FLAG2,max_iterMVI );

save([work_dir '\Mvec_TMVI_iter_.fld'],'-ascii','M')
write_MAG3D_TMI([work_dir '\TMVI_iter_.pre'],H,I,Dazm,obsx,obsy,obsz,pred_TMI,wd)

% Normalize magnetization vector
esus_MVI = sqrt( sum( M.^2 , 2 ) );
m = spdiags( 1./ esus_MVI ,0,mcell,mcell) * M;

% All zeros vector magnetization are put back to induced
indx = esus_MVI==0;

m(indx,:) = kron(ones(sum(indx),1),P');
fprintf('END AMI ITERATION # %i\n',ii);
end

