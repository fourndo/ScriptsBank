% Computes weighting matrix inside UBC - MAG3D inversion.
% Apply compactness by creating w.dat file after computing Ws, Wx, Wz
close all
clear all

addpath functions

meshfile = importdata('Mesh_10m.msh', ' ', 0);
meshfile=meshfile';
[dX,dY,dZ]=meshmaker(meshfile);

% data=importdata('Gravity_UC5m_Bcorr_RegRemov.obs', ' ', 3);
% target_misfit=str2num(data.textdata{3});
% d_obs=data.data(:,4);

% dos('cd C:\PVK_Projects\3796_Compactness\MAG3D');
% dos('maginv3d_newbounds inv_ctrl_round1.inp');

% [misfit]=get_misfit(d_obs);


nX=meshfile(1,1); %size(X,1);    %number of cell in X
nY=meshfile(2,1); %size(X,2);    %number of cell in Y
nZ=meshfile(3,1); %size(X,3);    %number of cell in Z



mcell=nX*nY*nZ;      %Number of model cells  

%Length scales
Lx=30;         
Ly=30;
Lz=30;

%Coefficient applied to WstWs in UBC
alphaS=1/Lx^2;


alphac=1.0e-0;  %Weight on compactness
alphax=1.0e-0;  %Coefficient for horizontal smoothness 
alphay=1.0e-0;
alphaz=1.0e-0;  %Coefficient for vertical smoothness

beta_in=2e+2;   %Trade-off parameter
delta=1e-11;     %Model disturbance 
lambda=1e-3;     %Scaling between smooth and compact


%Lp-norm
p=0;


% Ws=ones(i*j*k,1);
oo=1;

%% Loop over using initial model
while oo<=7%misfit>=(target_misfit/2)
wx=ones((nX-1)*nY*nZ,1)*alphax;
wy=ones(nX*(nY-1)*nZ,1)*alphay;
wz=ones(nX*nY*(nZ-1),1)*alphaz;    
  
% Import the file
    if oo==1
%         model=ones(i*j*k,1)*0.001;
        Wc=ones(mcell,1);
        beta(oo)=beta_in;
    else
    model = importdata('maginv3d_001.sus');
    
    beta(oo)=beta_in/4;

    gradm=comp_gradm(model, nX, nY, nZ, dX, dY, dZ);

%Small model weighting
    wc=1 ./ ( model .^ (2-p) + delta );
    
    phi_c=model'*(wc.*model);
    phi_m=model'*model;
    phi_xyz=model(model~=-1)' * gradm(model~=-1);
    
    alphac=phi_xyz / phi_c * lambda;

    Wc=alphac*wc;


    % phiWc(oo)=((Wc')*model.^2)*alphaS*Lx*Lz*Ly

    BB=sort(model(model~=-1));WctWc= 1.0 ./ ( BB .^ (2-p) + delta ) * alphac;
    WctWcm=(WctWc).*BB;
    figure(1);plot(BB,WctWcm,'r');
    xlim([0 max(model)])
    figure(2); hist(WctWcm(BB~=-1),50);
    figure(3); hist(model(model~=-1),30);
%     beta(oo)=beta_in%*(misfit(oo)/(misfit(1)));%/(2*oo);%         %Trade-off parameter
    end

    save('input_w.dat', 'Wc','-ascii')
    fid=fopen('input_w.dat','a');
    fprintf(fid,'\n');
    fclose(fid);

    save('input_w.dat', 'wx','-ascii','-append')
    fid=fopen('input_w.dat','a');
    fprintf(fid,'\n');
    fclose(fid);
    
    save('input_w.dat', 'wy','-ascii','-append')
    fid=fopen('input_w.dat','a');
    fprintf(fid,'\n');
    fclose(fid);
    

    save('input_w.dat', 'wz','-ascii','-append')
    fid=fopen('input_w.dat','a');
    fprintf(fid,'\n');
    fclose(fid);

    if oo==1
    fid=fopen('inv_ctrl.dat','w');
    fprintf(fid,'0 ! irest\n');
    fprintf(fid,'2 ! irest\n');
    fprintf(fid,'%12.8f 0.02   ! par tolc\n',beta(oo));
    fprintf(fid,'Obs_intrusive_2pc.obs  ! observations file\n');
    fprintf(fid,'maginv3d.mtx\n');
    fprintf(fid,'null  ! initial model\n');
    fprintf(fid,'null  ! reference model\n');
    fprintf(fid,'null  !active cell file\n');
    fprintf(fid,'null  ! lower, upper bounds\n');
    fprintf(fid,'%d %d %d  ! Le, Ln, Lz\n', Lx,Ly,Lz);
    fprintf(fid,'SMOOTH_MOD\n');
    fprintf(fid,'null  ! weighting file\n');
    fprintf(fid,'0');
    fclose(fid); 
    else
    %Update ctrl file
    fid=fopen('inv_ctrl.dat','w');
    fprintf(fid,'0 ! irest\n');
    fprintf(fid,'2 ! irest\n');
    fprintf(fid,'%12.8f 0.02   ! par tolc\n',beta(oo));
    fprintf(fid,'Obs_intrusive_2pc.obs  ! observations file\n');
    fprintf(fid,'maginv3d.mtx\n');
    fprintf(fid,'maginv3d_001.sus  ! initial model\n');
    fprintf(fid,'null  ! reference model\n');
    fprintf(fid,'null  !active cell file\n');
    fprintf(fid,'null  ! lower, upper bounds\n');
    fprintf(fid,'%d %d %d  ! Le, Ln, Lz\n',Lx,Ly,Lz);
    fprintf(fid,'SMOOTH_MOD\n');
    fprintf(fid,'input_w.dat  ! weighting file\n');
    fprintf(fid,'0');
    fclose(fid);
    end

oo=oo+1;

clear wx wy wz Wc WctWc WctWcm

% UBC Inversion
dos('maginv3d_newbounds inv_ctrl.dat');


end

% figure;plot(misfit)
% figure;plot(2:oo-1,phiWxyz(2:oo-1),2:oo-1,phiWc(2:oo-1),'r')