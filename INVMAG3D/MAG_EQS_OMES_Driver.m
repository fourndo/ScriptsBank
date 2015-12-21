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
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Equi_source\Tripple_Block_lined';
inpfile   = 'MAG_EQS_OMES.inp'; 

[meshfile,obsfile,topofile,mstart,mref,chi_target,alphas,beta,nlayer,FLAG1] = MAG3C_OMES_read_inp([work_dir '\' inpfile]);
FLAG2 = 'GRADm';

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

% Create nullcell
[nullcell,tcellID,ztopo_n] = topocheck(Xn,Yn,Zn,topo);
save([work_dir '\nullcell.dat'],'-ascii','nullcell');
% load([work_dir '\nullcell.dat']);
%% Load observation file (3C UBC-MAG format)
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);
% plot_mag3C(obsx,obsy,d,I,D,'Observed 3C-data')
% plot_TMI(obsx,obsy,d,d,wd,'Observed vs Predicted Magnitude');
% obsx = obsx(1:2:end);
% obsy = obsy(1:2:end);
% obsz = obsz(1:2:end);
% d = d(1:2:end);
% wd = wd(1:2:end);

% [XYZd_out] = Filter_xy([obsx obsy obsz d wd],20);

% obsx = XYZd_out(:,1);
% obsy = XYZd_out(:,2);
% obsz = XYZd_out(:,3);
% d = XYZd_out(:,4);
% wd = XYZd_out(:,5);

ndata = length(d);
Wd   = spdiags(1./wd,0,ndata,ndata);

nstn = length(obsx);

%% Get index of surface active cells
OMESid3D = MAG3C_get_OMES_id( reshape(nullcell,nz,nx,ny) , nlayer );
% index = obsy < (obsx + 9);
% OMESid3D = MAG3C_get_OMES_id_v2( obsx, obsy, xn, yn, zn, nullcell, nlayer );

% Convert from 3D to 2D mesh
[~,Ii,Jj] = ind2sub([nz nx ny],OMESid3D);

OMESid2D = sub2ind([nx ny],Ii,Jj);
mcell = 3*length(OMESid2D);

Tindex = zeros(1,nx*ny); Tindex(OMESid2D) = 1;
Tindex = kron(Tindex,ones(1,3));


% Get nodal discretization for octree levels
acelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,OMESid3D,ztopo_n,0);
% acelln = reshape( acelln , size(acelln,1) , 6 );

% Stretch the cells down
acelln(:,4) = acelln(:,4) - 1000;
dz = dz + 1000;
% acellnx = acelln;
% acellny = acelln;
% acellnz = acelln;

% Augmente the system by 3 for orthogonal components
% acelln = kron(acelln,ones(3,1));
% acelln = reshape( acelln ,size(acelln,1) , 1 , 6 );

% Create three orthogonal systems
% axc = (acellnx(:,2) + acellnx(:,5)) / 2;
% acellnx(:,2) = axc - dx(1)/20;
% acellnx(:,5) = axc + dx(1)/20;
% 
% ayc = (acellny(:,3) + acellny(:,6)) / 2;
% acellny(:,3) = ayc - dy(1)/20;
% acellny(:,6) = ayc + dy(1)/20;
% 
% azc = (acellnz(:,1) + acellnz(:,4)) / 2;
% acellnz(:,1) = azc + dz(1)/20;
% acellnz(:,4) = azc - dz(1)/20;
% 
% acelln(1:3:end,1,:) = acellnx;
% acelln(2:3:end,1,:) = acellny;
% acelln(3:3:end,1,:) = acellnz;
% % 
acelln = reshape( acelln ,size(acelln,1) , 1 , 6 );
%% Create model magnetization vectors
% azm = ones(3,1);
% azm(1) = 90;
% azm(2) = 0;
% azm(3) = 0;
% 
% dip = ones(3,1);
% dip(1) = 0;
% dip(2) = 0;
% dip(3) = 90;
% 
% mag_azmdip = [azm(:) dip(:)];
% mag_azmdip = kron( ones( mcell,1 ) , mag_azmdip);
% 
% mag_xyz = azmdip_2_xyz( mag_azmdip(:,1) , mag_azmdip(:,2) );

% save([work_dir '\OML.dat'],'-ascii','mag_azmdip');

% Create or load reference model
if ischar(mref)==1
    
    mref = load([work_dir '\' mref]);
    
else
    
%     mref = ones(mcell,1)*mref;
%     mref = [cosd(I)*cosd(D);(cosd(I)*sind(D));sind(I)]*1e-2;
    mref = zeros(mcell,1);

end

% Create or load reference model
if ischar(mstart)==1
    
    mstart = load([work_dir '\' mstart]);
    
else
    
    mstart = ones(mcell,1)*mstart;
    
end

%% Create model magnetization vectors
% Azimuth and dip of magnitization

% M = [spdiags(H * mag_xyz(:,1),0 , mcell , mcell);...
%      spdiags(H * mag_xyz(:,2),0 , mcell , mcell);...
%      spdiags(H * mag_xyz(:,3),0 , mcell , mcell)];


%% TMI forward projector
Ptmi = [spdiags(ones(nstn,1)* (cosd(I) * cosd(D)),0,nstn,nstn) ...
    spdiags(ones(nstn,1)* (cosd(I) * sind(D)),0,nstn,nstn) ...
    spdiags(ones(nstn,1)* sind(I),0,nstn,nstn)];

%% Compute depth weighting
% wr = get_wr(obsx, obsy, obsz, D, I, xn, yn, zn, nullcell, 'CENTER');
% save([work_dir '\wr.dat'],'-ascii','wr');

% wr = wr(OMESid3D);

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

    Tx(ii,:) = tx*H;
    Ty(ii,:) = ty*H;
    Tz(ii,:) = tz*H;
    
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

nactv = size(G,2)/3;

% clear Tx Ty Tz
wr = ones(size(G,2),1);
% wr(1:nactv) = sum(abs(G(:,(1):nactv)),1)./sum(abs(G(:,(1+2*nactv):3*nactv)),1);
% wr(1+nactv:2*nactv) = sum(abs(G(:,1+nactv:2*nactv)),1)./sum(abs(G(:,(1+2*nactv):3*nactv)),1);
% wr = 1./sum(abs(G),1)';%%ones(size(G,2),1);%
% wr = sqrt(wr ./ max(wr));
% wr(1+nactv:2*nactv) = sqrt(wr(1+nactv:2*nactv) ./ max(wr(1+nactv:2*nactv)));
% wr(1+2*nactv:3*nactv) = sqrt(wr(1+2*nactv:3*nactv) ./ max(wr(1+2*nactv:3*nactv)));

% wr_3D = zeros(nx*ny*nz,1); wr_3D(OMESid3D) = wr;
% 
% save([work_dir '\wr.dat'],'-ascii','wr_3D');

Wr = spdiags(wr,0,mcell,mcell);
IWr = spdiags(1./wr,0,mcell,mcell);



%% Create gradient matrices and corresponding volume vectors
nullcell = zeros(1,nx*ny);
nullcell(OMESid2D) = 1;
nullcell = (kron(nullcell,ones(1,nlayer)));

X = speye(nx*ny*nlayer);
X = X(nullcell==1,:);

[Ws, Wx, Wy, Wz, V, Vx, Vy, Vz] = get_GRAD_op3D_SQUARE(dx,dy,dz(1:nlayer),nullcell,X);

% [Ws, V ] = getWs3D(dx,dy,dz);
Ws = kron(speye(3), Ws);
Wx = kron(speye(3), Wx);
Wy = kron(speye(3), Wy);
Wz = kron(speye(3), sparse(size(Wz,1),size(Wz,1)));
Vx = kron(speye(3), Vx);
Vy = kron(speye(3), Vy);
Vz = kron(speye(3), Vz);
V = kron( speye(3), V );


%% Apply depth and data weighting on sensitivity
% sens = mean((abs(G*V)),1)'/ndata;
% save([work_dir '\sens.dat'],'-ascii','sens');
G   = Wd * G;
d = Wd * d;

%% Inversion   
delta=1e-10;     %HARDWIRED: Small term in compact function

target = chi_target * ndata;     % Target misifit

% Compute total objective function
comp_phi = @(m,phim,b) sum( ( G * m - d ).^2 ) + (m)' * b * phim * (m);

% Wx = Wx *Wr;
% Wy = Wy *Wr;
% Wz = Wz *Wr;

counter = 1;               

                
% Initialize inversion
invmod      = mstart;       % Initial model       

phi_init    = sum((G * invmod - d).^2);   % Initial misfit
phi_d       = phi_init;
phi_m       = [];         
lp_count = 0;
LP = [2 0 0 1];
% alphas(4) = 0;
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
lowBvec = ones(mcell,1) * lowB;
uppBvec = ones(mcell,1) * uppB;
objfunc = @(m,phim,b) sum( ( G * m - d ).^2 ) + (m)' * b * phim * (m);
if nlayer ==1
    
    aVRWs = alphas(1) * ( V * Ws );
    aVRWx = alphas(2) * ( Vx * Wx );
    aVRWy = alphas(3) * ( Vy * Wy );
    aVRWz = alphas(4) * ( Vz * Wz );
    
else
    
    aVRWs = alphas(1) * ( V * Ws );
    aVRWx = alphas(2) * ( Vx * Wx );
    aVRWy = alphas(3) * ( Vy * Wy );
    aVRWz = alphas(4) * ( Vz * Wz );
    
end

MOF = aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWy'*aVRWy + aVRWz'*aVRWz;
ncg = 0;
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

%         phi = norm(G*invmod - d).^2 +...
%             invmod' * beta * MOF * invmod;

phi(count) = objfunc(invmod,MOF,beta(count));
        phi_m(count) = (invmod)'*(MOF)*(invmod);
        delta(count) = 1;
        L2 = [2 2 2 1];
        [MOF,aVRWs,aVRWx,aVRWy] =  get_lp_MOF_2D(invmod,mref,phi_m(:,end),ones(mcell,1),ones(mcell,1),V,Ws,Vx,Wx,Vy,Wy,wr,alphas,LP,FLAG2,FLAG1,0,delta(count));

        
    else
        lp_count = lp_count+1;
                        
        if delta(end) > 1e-3
            
            delta(count) = delta(count-1)/10;

        else

            delta(count) = delta(count-1);

        end


         fprintf('\n# # LP-LQ ITER# #\n');
        [MOF,aVRWs,aVRWx,aVRWy] = get_lp_MOF_2D(invmod,mref,phi_m(:,end),ones(mcell,1),ones(mcell,1),V,Ws,Vx,Wx,Vy,Wy,wr,alphas,LP,FLAG2,FLAG1,1,delta(count));

    end

    

    %% Pre-conditionner
    diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
    PreC     = Pac * spdiags(1./diagA(:),0,mcell,mcell);

    m_in = invmod;
    
    [invmod, iter, Pac] = GN_PCG_solver( G, invmod, zeros(mcell,1), nullcell, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );

    %% Step length, line search                
    ncg = ncg+iter; % Record the number of CG iterations  

    %% Save iteration and continue
    clear A

    phi(count) = comp_phi(invmod,MOF,beta(count));

    phi_d(count) = sum(( G*invmod - d ).^2);

    % Cool beta
    if phi_d(count) < target*2 && count~=1

      beta(count+1) = 0.5 * beta(count);

    else

      beta(count+1) = 0.5 * beta(count);

    end

    phi(count) = objfunc(invmod,MOF,beta(count));
    phi_m(count) = (invmod)'*(MOF)*(invmod);
    fprintf(fid,' \t %i \t %8.5e ',count,beta(count));
    fprintf('Iteration: \t %i  \nBeta: \t %8.5e \n',count,beta(count));
    fprintf(fid,' \t %8.5e ',phi_d(count));
    fprintf('phid:\t %8.5e\n',phi_d(count))
    fprintf(fid,' \t %8.5e ',sum((aVRWs*invmod).^2));
    fprintf(fid,' \t %8.5e ',sum((aVRWx*invmod).^2));
    fprintf(fid,' \t %8.5e ',sum((aVRWy*invmod).^2));
    fprintf(fid,' \t %8.5e ',sum((aVRWz*invmod).^2));
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
m2D = invmod;

mcell = length(OMESid2D);

kx = zeros(nx*ny*nz,1); kx(OMESid3D) = m2D(1:mcell);
ky = zeros(nx*ny*nz,1); ky(OMESid3D) = m2D(1+mcell:2*mcell);
kz = zeros(nx*ny*nz,1); kz(OMESid3D) = m2D(1+2*mcell:3*mcell);

eqs_3D = zeros(nx*ny*nz,1); eqs_3D(OMESid3D) = sqrt(m2D(1:mcell).^2 + m2D(1+mcell:2*mcell).^2 + m2D(1+2*mcell:3*mcell).^2);
save([work_dir '\EQS_3D.sus'],'-ascii','eqs_3D');
% save([work_dir '\OMES_kx.sus'],'-ascii','kx');
% save([work_dir '\OMES_ky.sus'],'-ascii','ky');
% save([work_dir '\OMES_kz.sus'],'-ascii','kz');

% data_3C = G * invmod;

% bx = G(:,1:mcell) * invmod(1:mcell);
% by = G(:,mcell+1:2*mcell) * invmod(mcell+1:2*mcell); 
% bz = G(:,2*mcell+1:3*mcell) * invmod(2*mcell+1:3*mcell); 

bx = Tx*(m2D); 
by = Ty*(m2D);
bz = Tz*(m2D);

write_MAG3D_TMI([work_dir '\OMES_TMI.pre'],H,I,Dazm,obsx,obsy,obsz,Ptmi * [bx;by;bz],wd);

write_MAG3D_TMI([work_dir '\OMES_lBl.pre'],H,I,Dazm,obsx,obsy,obsz,sqrt(bx.^2 + by.^2 + bz.^2),wd);

write_MAG3D_3C([work_dir '\OMES_3C.pre'],H,I,Dazm,obsx,obsy,obsz,bx,by,bz,ones(ndata,1),ones(ndata,1),ones(ndata,1));

plot_mag3C(obsx,obsy,[bx;by;bz],I,D,'Observed 3C-Data data')

fprintf(fid,'End of OMES inversion. Number of iterations: %i\n',count);
%             fprintf('Final data misfit: %8.3e. Final l1-model error: %8.3e\n\n',phi_d(count),norm(m-model_out,1))
fclose(fid);


lMl = [kx ky kz];
save([work_dir '\OMES_kvec.fld'],'-ascii','lMl');

%% Forward model RTP
% Rz = @(x)   [cosd(x) -sind(x) 0;
%             sind(x) cosd(x) 0;
%             0 0 1];
% Rx = @(x)   [1 0 0;
%             0 cosd(x) -sind(x);
%             0 sind(x) cosd(x)]; 
%                 
% Ry = @(x)   [cosd(x) 0 sind(x);
%             0 1 0;
%             -sind(x) 0 cosd(x)];
% 
% RR = kron(speye(nx*ny),Rx(90) *Rz(0));
% d_RTP = Tz * (RR * model_out);
%%
% Get nodal discretization for octree levels
% mOMES = size(OMESid2D,1);
% acelln = MAG3C_RTC_OctNodes(Xn,Yn,Zn,OMESid2D,ztopo_n,0);
% acelln = reshape( acelln , mOMES , 6 );
% 
% % Stretch the cells down
% acelln(:,1) = acelln(:,1) - 1000;
% dz = dz + 1000;
% 
% acelln = reshape( acelln , mcell , 1 , 6 );


% rtptmi = [(cosd(90) * cosd(90)) (cosd(90) * sind(90)) sind(90)];

% model_out = sqrt(m2D(1:3:end).^2 + m2D(2:3:end).^2 + m2D(3:3:end).^2 );

% m_3D_out = zeros(nz,nx,ny);
% m_3D_out(1,:,:) = reshape(model_out,nx,ny);m_3D_out = m_3D_out(:);
% save([work_dir '\OMES_Ampl.sus'],'-ascii','m_3D_out');

% d_RTP = zeros(ndata,1);
% progress = -1;
% tic  
% for ii = 1:ndata
% 
%     
%     % compute kernel for active cells
%     [tx,ty,tz] = MAG3C_T(obsx(ii),obsy(ii),obsz(ii),acelln);
% 
%     d_RTP(ii) = tz*M*model_out;
% 
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
