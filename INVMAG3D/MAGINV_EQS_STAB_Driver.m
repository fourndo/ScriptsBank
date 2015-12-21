% ORTHOGONALLY MAGNETIZED EQUIVALENT SOURCE 
% Written by: Dominique Fournier 
% Last update: July 30th, 2014 
% 
% The code generates an equivalent source layer along topography. The layer
% can be use to forward model components and predict magnetization
% direction


clear all
close all

uppB = 10;
lowB = -10;

addpath ..\FUNC_LIB\;
% Project folders
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Equi_source\SingleBlock_Li\Corner_off';
inpfile   = 'MAG_EQS_OMES.inp'; 

[meshfile,obsfile,topofile,mstart,mref,chi_target,alphas,beta,nlayer,FLAG1] = MAG3C_OMES_read_inp([work_dir '\' inpfile]);

% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

% Write logfile
fid = fopen([work_dir '\MAG3Csen.log'],'w');
fprintf(fid,'MAG3Csen\n');
fprintf(fid,'Generates sparse matrices for magnetostatic forward modeling: Tx, Ty, Tz\n');
fprintf(fid,'Topographic model: nullcell.dat\n');
fprintf(fid,'DISTANCE | DEPTH weighting: wr.dat\n\n');
fprintf(fid,'Written by: Dominique Fournier\n');
fprintf(fid,'Last update: July 14th, 2014\n\n');
fprintf(fid,'INPUT FILES: \n');
fprintf(fid,'Mesh: \t\t\t %s \n',meshfile);
fprintf(fid,'Obsfile: \t\t\t %s \n',obsfile);
fprintf(fid,'Topography: \t\t\t %s \n',topofile);
fclose(fid);

%% 3D nodal location
[Zn,Xn,Yn] = ndgrid(zn,xn,yn);

if isempty(topofile)==1
    
    nxn = length(xn);
    nyn = length(yn);
    topo = [reshape(Xn(1,:,:),nxn*nyn,1) reshape(Yn(1,:,:),nxn*nyn,1) reshape(Zn(1,:,:),nxn*nyn,1)+1e-8];
    
else
    % Load topo
    topo = read_UBC_topo([work_dir '\' topofile]);
    

end

%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, data, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');

ndat_orig = length(data);
ndata = length(data);
Wd   = spdiags(1./wd,0,ndata,ndata);

nstn = length(obsx);


%% Manually add observations - TO BE AUTOMATED
% xx = -31.5:1:31.5;
% [YY,XX] = ndgrid(xx,xx);
% 
% index = find((YY(:) < XX(:)+10) & (YY(:) > XX(:)+0));
% padloc = [XX(index) YY(index) zeros(length(index),1)];
% 
% % Select edge of data
% % index_A = find((YY(:) < XX(:)+10) & (YY(:) > XX(:)+8));
% % index_B = find((YY(:) > XX(:)+30) & (YY(:) < XX(:)+32));
% % 
% lim_data = zeros(length(index),1);
% 
% for ii = 1 : length(index)
%     
%     lim_data(ii) = data(obsx==XX(index(ii)) & obsy==YY(index(ii)));
% 
% end
% % 
% % lim_data = [lim_data;zeros(length(index_B),1)];
% % 
% % F = scatteredInterpolant([XX(index_A);XX(index_B)],[YY(index_A);YY(index_B)],lim_data);
% % pad_data = F(XX(index),YY(index));
% 
% % pad_data = data(index);
% 
% obsx = [obsx;YY(index)-9.5];
% obsy = [obsy;XX(index)+9.5];
% data = [data;lim_data];
% obsz = [obsz;zeros(length(index),1)];
% wd = [wd;ones(length(index),1)];
% 
% index = obsx > -32 & obsy < 32;
% 
% data = data(index);
% obsx = obsx(index);
% obsy = obsy(index);
% obsz = obsz(index);
% wd = wd(index);
% 
% figure;scatter(obsx,obsy,10,data)
ndata = length(data);
Wd   = spdiags(1./wd,0,ndata,ndata);

nstn = length(obsx);


%% Extract active cells
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
% save([work_dir '\nullcell.dat'],'-ascii','nullcell');

% Get index of surface active cells
% OMESid = MAG3C_get_OMES_id( reshape(nullcell,nz,nx,ny) , nlayer );
% index = obsy < (obsx + 9);
OMESid = MAG3C_get_OMES_id_v2( obsx, obsy, xn, yn, zn, nullcell, nlayer);
[~,Ii,Jj] = ind2sub([nz nx ny],OMESid);

OMESid2D = sub2ind([nx ny],Ii,Jj);
Tindex = zeros(1,nx*ny); Tindex(OMESid2D) = 1;

mOMES = size(OMESid,1);

mcell = mOMES;

% Get nodal discretization for octree levels
acelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,OMESid,ztopo_n,0);
acelln = reshape( acelln , mOMES , 6 );

% Stretch the cells down
% acelln(:,4) = acelln(:,4)-500 ;
% dz = dz+500 ;

acelln = reshape( acelln , mcell , 1 , 6 );

%% Create model magnetization vectors
% azm = ones(3,1);
azm(1) = Dazm;
% azm(2) = 0;
% azm(3) = 0;

% dip = ones(3,1);
dip(1) = I;
% dip(2) = 90;
% dip(3) = 90;

% mag_azmdip = [azm(:) dip(:)];
mag_azmdip = kron( ones( mOMES,1 ) , [Dazm I]);

mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );




% save([work_dir '\OML.dat'],'-ascii','mag_azmdip');

% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir '\' mref]);
    
else
    
    mref = ones(mcell,1)*mref;
    
end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir '\' mstart]);
    
else
    
    mstart = ones(mcell,1)*mstart;
    
end

%% Create model magnetization vectors
% Azimuth and dip of magnitization

M = [spdiags(H * mag_xyz(:,1),0 , mcell , mcell);...
     spdiags(H * mag_xyz(:,2),0 , mcell , mcell);...
     spdiags(H * mag_xyz(:,3),0 , mcell , mcell)];


%% TMI forward projector
Ptmi = [spdiags(ones(nstn,1)* (cosd(I) * cosd(D)),0,nstn,nstn) ...
    spdiags(ones(nstn,1)* (cosd(I) * sind(D)),0,nstn,nstn) ...
    spdiags(ones(nstn,1)* sind(I),0,nstn,nstn)];



%% Generate sensitivities matrix generated by FMAG3C
% Pre-allocate to store fields
Tx = zeros(ndata,mcell);
Ty = zeros(ndata,mcell);
Tz = zeros(ndata,mcell);

progress = -1;
tic  
for ii = 1:ndata

    
    % compute kernel for active cells
    [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),acelln);

    Tx(ii,:) = tx*M;
    Ty(ii,:) = ty*M;
    Tz(ii,:) = tz*M;
    
    d_iter = floor(ii/ndata*20);
    if  d_iter > progress

        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
        progress = d_iter;

    end
            
end

% save([work_dir '\OMESTx'],'Tx');
% save([work_dir '\OMESTy'],'Ty');
% save([work_dir '\OMESTz'],'Tz');
% 
% load([work_dir '\OMESTx']);
% load([work_dir '\OMESTy']);
% load([work_dir '\OMESTz']);

G = Ptmi*[Tx;Ty;Tz];

wr = ones(size(G,2),1);%sum(abs(G),1)';
wr = sqrt(wr ./ max(wr));

% wr_3D = zeros(nx*ny*nz,1); wr_3D(OMESid) = wr;

% save([work_dir '\wr.dat'],'-ascii','wr_3D');

Wr = spdiags(wr,0,mcell,mcell);
IWr = spdiags(1./wr,0,mcell,mcell);
% clear Tx Ty Tz

% G = G;

%% Compensate each uncertainties for edges effect

% avgG = mean(G,1);
% scale = zeros(ndata,1);
% xc = mean(acelln(:,1,[2 5]),3);
% yc = mean(acelln(:,1,[3 6]),3);
% for ii = 1 : length(obsx)
%     
%     
%     % Find closest cells
%     R = sqrt( (obsx(ii) - xc).^2 + (obsy(ii) - yc).^2);
%     
%     index = find(R==min(R));
%     
%     scale(ii) = abs(avgG(index ));
%     
% end
% 
% scale = scale / max(scale);
% 
% figure; scatter(obsx,obsy,10,scale);
% colorbar
% 
% wd = wd.*scale;
% Wd   = spdiags(1./wd,0,ndata,ndata);

%% Create gradient matrices and corresponding volume vectors
% [Wx, Wy, Wz, Vx, Vy, Vz] = get_GRAD_op3D_v4(dx,dy,dz(1:nlayer)',ones(mcell,1));
% [Ws, v ] = getWs3D(dx,dy,dz(1:nlayer)');
nullcell = zeros(1,nx*ny*nlayer);
nullcell(OMESid) = 1;

X = speye(nx*ny*nlayer);
X = X(OMESid,:);

[Ws, Wx, Wy, Wz, V, Vx, Vy, Vz] = get_GRAD_op3D_SQUARE(dx,dy,dz(1:nlayer),nullcell,X);

% Boost smoothing
% booster = Tindex(1:end-1)==1 & Tindex(2:end)==0;
% booster = [booster 0];
% 
% Wx(booster==1,:) = Wx(booster==1,:)*100;
% Wy(booster==1,:) = Wy(booster==1,:)*100;

% Ws = Ws(Tindex==1,Tindex==1);
% Wx = Wx(Tindex==1,Tindex==1);
% Wy = Wy(Tindex==1,Tindex==1);
% Wz = Wz(Tindex==1,Tindex==1);
% V = V(Tindex==1,Tindex==1);
% Vx = Vx(Tindex==1,Tindex==1);
% Vy = Vy(Tindex==1,Tindex==1);
% Vz = Vz(Tindex==1,Tindex==1);

%% Apply depth and data weighting on sensitivity
% sens = mean((abs(G*V)),1)'/ndata;
% save([work_dir '\sens.dat'],'-ascii','sens');
G = G * IWr;
G   = Wd * G;
d = Wd * data;

%% Inversion   
delta=1e-10;     %HARDWIRED: Small term in compact function

target = chi_target * ndata;     % Target misifit

% Compute total objective function
comp_phi = @(m,phim,b) sum( ( G * m - d ).^2 ) + (m)' * b * phim * (m);

% Wx = Wx *Wr;
% Wy = Wy *Wr;
% Wz = Wz *Wr;
lowBvec = ones(mcell,1) * lowB;
uppBvec = ones(mcell,1) * uppB;
counter = 1;               
ncg = 0;
                
% Initialize inversion
invmod      = Wr*mstart;       % Initial model       

phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phi_init;
phi_m       = [];         

count=1;
countbar = 0;

% Message prompt
logfile = [work_dir '\Log_OMS.log'];
fid = fopen(logfile,'w');
fprintf(fid,'Starting OMES inversion\n');
fprintf(fid,'Starting misfit %e\n',phi_init);
fprintf(fid,'Target misfit %e\n',target);
fprintf(fid,'Iteration:\t\tBeta\t\tphid\t\tphis\t\t ');
fprintf(fid,'phix\t\tphiy\t\tphiz\t\tphim\t\tphi\t\t ');
fprintf(fid,'#cut cells \t # CG Iter\n');


% Initiate active cell
lowb = zeros(mcell,1); % Logical for lower boundded cells
uppb = zeros(mcell,1); % Logical for upper boundded cells
Pac = spdiags((lowb==0).*(uppb==0),0,mcell,mcell);

WxtWx = ( Vx * Wx )' * ( Vx * Wx ) ;
WytWy = ( Vy * Wy )' * ( Vy * Wy ) ;
WztWz = ( Vz * Wz )' * ( Vz * Wz ) ;

WstWs = ( V * Ws )' * ( V * Ws ) ;

if nlayer ==1
    
    aVRWs = alphas(1) * ( V * Ws );
    aVRWx = alphas(2) * ( Vx * Wx );
    aVRWy = alphas(3) * ( Vy * Wy );
    aVRWz = sparse ( mcell , mcell );
    
else
    
    aVRWs = alphas(1) * ( V * Ws );
    aVRWx = alphas(2) * ( Vx * Wx );
    aVRWy = alphas(3) * ( Vy * Wy );
    aVRWz = alphas(4) * ( Vz * Wz );
    
end

if nlayer ==1
    
    MOF = alphas(1)*WstWs + alphas(2)*WxtWx + alphas(3)*WytWy;
    
else
    
    MOF = alphas(1)*WstWs + alphas(2)*WxtWx + alphas(3)*WytWy + alphas(4)*WztWz;
    
end

while phi_d(end) > target 


    %First iteration uses the usual smallness and smoothness
    %regularization. Choose a beta on the ratio of the traces 
    % of the elements of objective function
    if count==1                 

        if isempty(beta)==1
            
            m_temp = randn(mcell,1);
%             beta = full( sum(sum(G.^2,1)) / sum(diag(MOF,0)) * 1e+1 );

            beta = sum((G*m_temp).^2) / (m_temp'*MOF*m_temp) * 1e+2;
        end
%                     beta = 1e+5;

        lambda = beta ;
        phi = norm(G*invmod - d).^2 +...
            invmod' * lambda * MOF * invmod;

    else

        lambda = beta(count);

    end


    %% Pre-conditionner
    diagA = sum(G.^2,1) + lambda*spdiags(MOF,0)';
    PreC     = Pac * spdiags(1./diagA(:),0,mcell,mcell);

    m_in = invmod;
    
    [invmod, iter, Pac] = GN_PCG_solver( G, invmod, mref, nullcell, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );

    
    %% Pre-conditionner
%     diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
%     PreC     = spdiags(1./diagA(:),0,mcell,mcell);

    %% PCG steps
%     A= [ G   ;...
%     sqrt( beta(count) ) * ( aVRWs ) ;...
%     sqrt( beta(count) ) * ( aVRWx ) ;...
%     sqrt( beta(count) ) * ( aVRWy ) ];
%         
% 
%    
%     b = [ d ; ...
%         zeros(size(aVRWs,1),1) ;...
%         zeros(size(aVRWx,1),1) ;...
%         zeros(size(aVRWy,1),1) ];
%     [invmod,~,iter] = PCGLSQ( invmod, A , b, PreC, Pac );
    %% Step length, line search                
    ncg = ncg+iter; % Record the number of CG iterations

    

    %% Save iteration and continue
    clear A

    phi(count) = comp_phi(invmod,MOF,lambda);

    phi_d(count) = sum(( G*invmod - d ).^2);

    % Cool beta
    if phi_d(count) < target*2 && count~=1

      beta(count+1) = 0.5 * beta(count);

    else

      beta(count+1) = 0.5 * beta(count);

    end

    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf('Iteration: \t %i  \nBeta: \t %8.5e \n',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf('phid:\t %8.5e\n',phi_d(count))
    fprintf(fid,' \t %8.5e ',invmod'*alphas(1)*WstWs*invmod);
    fprintf(fid,' \t %8.5e ',invmod'*alphas(2)*WxtWx*invmod);
    fprintf(fid,' \t %8.5e ',invmod'*alphas(3)*WytWy*invmod);
    fprintf(fid,' \t %8.5e ',invmod'*alphas(4)*WztWz*invmod);
    fprintf(fid,' \t %8.5e ',invmod'*MOF*invmod);
    fprintf(fid,' \t %8.5e ',phi(count));
    fprintf(fid,' \t\t %i ',sum(lowb));
    fprintf('Number of Inactive cells: %i\n',sum(lowb));
    fprintf(fid,' \t\t %i\n',ncg);
    fprintf('Number of CGS iterations: %i\n',ncg);
   % Output interation result

%                 save([work_dir '\OMES_' 'p' num2str(pvec(pp)) 'q' num2str(qvec(qq)) 'l' num2str(lvec(ll)) '_iter_' num2str(count) '.sus'],'-ascii','model_out')
%                 write_MAG3D_TMI([work_dir '\TMI_iter_' num2str(count) '.pre'],H,I,Dazm,obsx,obsy,obsz,(G*invmod).*wd,wd)

    count=count+1;
    countbar = 0;  
    
end
%%
count=count-1;  
invmod = IWr*invmod;

eqs_3D = zeros(nx*ny*nz,1); eqs_3D(OMESid) = invmod;

save([work_dir '\EQS_3D.sus'],'-ascii','eqs_3D');

S = spdiags(ones(ndat_orig,1),0,ndat_orig,ndata);
d_pre = S * G * invmod;
wd = S * wd;
obsx = S * obsx;
obsy = S * obsy;
obsz = S * obsz;

plot_TMI(obsx,obsy,S*data,d_pre.*wd,ones(ndat_orig,1),'Observed vs Predicted Magnitude');

bx = S * Tx*invmod;%data_3C(1:ndata); 
by = S * Ty*invmod;%data_3C(ndata+1:2*ndata); 
bz = S * Tz*invmod;%data_3C(2*ndata+1:3*ndata); 



write_MAG3D_TMI([work_dir '\OMES_TMI.pre'],H,I,Dazm,obsx,obsy,obsz,d_pre,wd);

write_MAG3D_TMI([work_dir '\OMES_lBl.pre'],H,I,Dazm,obsx,obsy,obsz,sqrt(bx.^2 + by.^2 + bz.^2),wd);

write_MAG3D_3C([work_dir '\OMES_3C.pre'],H,I,Dazm,obsx,obsy,obsz,bx,by,bz,ones(ndat_orig,1),ones(ndat_orig,1),ones(ndat_orig,1));

plot_mag3C(obsx,obsy,[bx;by;bz],I,D,'Observed 3C-Data data')

fprintf(fid,'End of OMES inversion. Number of iterations: %i\n',count);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
fclose(fid);


% lMl = [kx ky kz];
% save([work_dir '\OMES_kvec.fld'],'-ascii','lMl');


%% Compute RTP
% mag_azmdip = kron( ones( mOMES,1 ) , [0 90]);
% 
% mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );
% 
% 
% % M = [spdiags(H * mag_xyz(:,1),0 , mcell , mcell);...
% %      spdiags(H * mag_xyz(:,2),0 , mcell , mcell);...
% %      spdiags(H * mag_xyz(:,3),0 , mcell , mcell)];
% M = [spalloc(mcell,mcell,0);...
%     spalloc(mcell,mcell,0);...
%     spdiags(H * ones(mcell,1),0,mcell,mcell)];
% % rtptmi = [(cosd(90) * cosd(90)) (cosd(90) * sind(90)) sind(90)];
% 
% % model_out = sqrt(m2D(1:3:end).^2 + m2D(2:3:end).^2 + m2D(3:3:end).^2 );
% 
% % rtptmi = [(cosd(90) * cosd(90)) (cosd(90) * sind(90)) sind(90)];
% 
% % Pre-allocate to store fields
% d_RTP = zeros(ndata,1);
% 
% progress = -1;
% tic  
% for ii = 1:ndata
% 
%     
%     % compute kernel for active cells
%     [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),acelln);
% 
%     d_RTP(ii) = tz*M*invmod;
%     
%     d_iter = floor(ii/ndata*20);
%     if  d_iter > progress
% 
%         fprintf('Computed %i pct of data in %8.5f sec\n',d_iter*5,toc)
%         progress = d_iter;
% 
%     end
%             
% end
% 
% write_MAG3D_TMI([work_dir '\' obsfile(1:end-4) '_RTP.pre'],H,I,Dazm,...
%     obsx,obsy,obsz,d_RTP,ones(ndata,1));
% 
% plot_TMI(obsx,obsy,d,d_RTP,ones(ndata,1),'Observed vs Predicted Magnitude');
