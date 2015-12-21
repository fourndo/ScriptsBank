% MAG_lp_Driver
% Author: D Fournier
% Last Update: April 21st, 2013
%
% This program invert for "compact" model using a speudo lp-norm on the
% model and derivative. It is an adapted algorithm based on the classic 
% paper of Last & Kubic. The central solver is the UBC - MAGINV3D code. 
%
% Apply compactness by creating a w.dat file, after computing Ws, Wx, Wz
% values. 
%
% INPUT
% Directory for the location of the following input files:
% maginv3d_001.sus: First outputed model from a unconstrained inversion.
% mesh_file.msh: Mesh file in UBC format
% topofile.topo: Topo file used for the unconstrained inversion.
% input.inp: Input file used for the unconstrained inversion.
% maginv3d.log: Log file from unconstrained inversion
%
% (for "continuous" mode only) 
% maginv3d.mtx: sensitivity matrix computed from magsen3d
%
% OUTPUT
% (for "step" mode)
% w.dat: Weighting matrix computed
% input.inp: Input file for the maginv3d
%
% (for "continuous" mode)
% magmodel.sus: Inverted mag models

addpath functions

close all
clear all

%% USER INPUTS
set_path = 'C:\Projects\3796_Compactness\MAG3D\Current\Input_test'   ;

dos (['copy maginv3d_001.sus' ' ' set_path]);
dos (['copy maginv3d.log' ' ' set_path]);

cd (set_path);

file_list=ls;

argin = 'continuous';


% Cycle through all the files to look for the inputs
input = [];
mesh = [];
logfile = [];
topo = [];
topo_model = [];
topo_check = [];

    for ii = 1:size(file_list,1)-2;

        look_at = strtrim(file_list(ii+2,:));

            if strcmp(look_at(end-5:end),'3d.log')==1

                logfile =look_at;

            elseif strcmp(look_at(end-2:end),'msh')==1

                mesh = look_at;

            elseif strcmp(look_at(end-2:end),'inp')==1

                input = look_at;

            elseif strcmp(look_at(end-3:end),'topo')==1

                topo =  look_at;

            elseif strcmp(look_at,'topo_model.txt')==1
            
                topo_check =  look_at;
            end
    end
    
    
    if isempty(topo_check)==1
        
    % Run a topo check if not done already
    fprintf('Computing Topocheck - might take few minutes')
    dos (['topocheck ' mesh ' ' topo])
    end
    
    load topo_model.txt
    
    

%% Read information from log file
    fid=fopen(logfile,'rt');
    
    misfit=[];
    max_num_lines = 30000;
    % Go through the log file and extract data and the last achieved misfit
    for ii=1:max_num_lines         	
    line=fgets(fid); %gets next line 

        if line==-1
            fprintf('File ended at line %i\n',ii);
            fprintf('Did not find the information needed - review log file\n')
            break
        end
        
        if length(strtrim(line))>=length('# of data:')
            description = strtrim(line);
            if strcmp(description(1:10),'# of data:')==1
                ndata = str2num(description(11:end));
            end
        end
        
        if length(strtrim(line))>=length('multiplier:')
            description = strtrim(line);
            if strcmp(description(1:11),'multiplier:')==1
                beta_in = str2num(description(12:end));
            end
        end
        
        % Extract alpha values
        if length(strtrim(line))>=length('Le, Ln, Lz:')
            description = strtrim(line);
            if strcmp(description(1:11),'Le, Ln, Lz:')==1
                lengths = str2num(description(12:end));
            end
        end
        
        if length(strtrim(line))>=length('data misfit:')
            description = strtrim(line);
            if strcmp(description(1:12),'data misfit:')==1
                misfit = str2num(description(13:end));
            end
        end
        
        if isempty(misfit)==0
            break
        end
        
    end

%% Extract mesh from UBC mesh file

addpath functions
meshfile = importdata(mesh, ' ', 0);

cd ..
meshfile=meshfile';
[dX,dY,dZ]=meshmaker(meshfile);

nX=meshfile(1,1); %size(X,1);    %number of cell in X
nY=meshfile(2,1); %size(X,2);    %number of cell in Y
nZ=meshfile(3,1); %size(X,3);    %number of cell in Z



mcell=nX*nY*nZ;      %Number of model cells  

% %Length scales
Lx=lengths(1);         
Ly=lengths(2);
Lz=lengths(3);

% %Coefficient applied to WstWs in UBC
alphaS=1/Lx^2;


alphac=alphaS;  %Weight on compactness
alphax=1.0e-0;  %Coefficient for horizontal smoothness 
alphay=1.0e-0;
alphaz=1.0e-0;  %Coefficient for vertical smoothness

% beta_in=2.5e+1;   %Trade-off parameter
delta=1e-11;     %Model disturbance 
lambda=1e-0;     %Scaling between smooth and compact
p = 0.1;
q = 2;


new_dir=(['Inv_lp' num2str(p*100) '_lq' num2str(q)]);
    
dos (['mkdir ' set_path '\' new_dir]);
%Weight on smoothness


% Ws=ones(i*j*k,1);
oo=1;
addpath functions
addpath Input_test

% Create derivative matrices
mnull = []; % No null values for now
% tic
% Wx = getWx_3D(mcell,dX,dY,dZ,mnull);
% Wy = getWy_3D(mcell,dX,dY,dZ,mnull);
% Wz = getWz_3D(mcell,dX,dY,dZ,mnull);
% toc
% WxtWx = Wx'*Wx;
% WytWy = Wy'*Wy;
% WztWz = Wz'*Wz;

wx=ones((nX-1)*nY*nZ,1)*alphax;
wy=ones(nX*(nY-1)*nZ,1)*alphay;
wz=ones(nX*nY*(nZ-1),1)*alphaz;  

Ws = getWs_3D(mcell,dX,dY,dZ,topo_model);
% clear Wx Wy Wz


%% Loop over using initial model
while oo <= 7 && misfit >= ndata
  
  
% Import the file

    model = load('maginv3d_001.sus');
    
    beta(oo)=beta_in/(1.5^(oo-1));

%     phi_xyz = model' * (alphax * WxtWx + alphay * WytWy + alphaz * WztWz) * model;
    Wxm=comp_gradxm(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
    Wym=comp_gradym(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
    Wzm=comp_gradzm(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
%     save('Wzm','-ascii','Wzm')
%     save('Wym','-ascii','Wym')
%     save('Wxm','-ascii','Wxm')
    phi_xyz = alphaz * Wzm' * Wzm + alphay * Wym' * Wym + alphax * Wxm' * Wxm;
    
    %Small model weighting
    wc=1 ./ ( (model.*topo_model) .^ (2-p) + delta );
    
    phi_c=(Ws.*model)'*(wc.*(Ws.*model));
%     phi_m=model'*model;
%     phi_xyz=model(model~=-1)' * gradm(model~=-1);
    
    alphac=phi_xyz/(phi_c*alphaS);

    Wc=alphac*wc;


    % phiWc(oo)=((Wc')*model.^2)*alphaS*Lx*Lz*Ly

    BB=sort(model(model~=-1));WctWc= 1.0 ./ ( BB .^ (2-p) + delta ) * alphac;
    WctWcm=(WctWc).*BB;
    figure(1);plot(BB,WctWcm,'r');
    xlim([0 max(model)])
    figure(2); hist(WctWcm(BB~=-1),50);
    figure(3); hist(model(model~=-1),30);
%     beta(oo)=beta_in%*(misfit(oo)/(misfit(1)));%/(2*oo);%         %Trade-off parameter

cd ([set_path]);
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
    
% Create input file

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

oo=oo+1;

% clear wx wy wz Wc WctWc WctWcm
switch argin 
    case 'continuous'
        
        % Call maginv3d_newbounds
       dos('maginv3d_newbounds inv_ctrl.dat');
       
       dos (['copy maginv3d_001.sus' ' ' set_path '\' new_dir])
       
       dos (['rename ' set_path '\' new_dir '\maginv3d_001.sus maginv3d_00' num2str(oo) '.sus '])
       
       cd ..
   case 'step-by-step'
       
    otherwise
        
        fprintf('Input argument not recognize\n')
        fprintf('Should be "step-by-step" or "continuous"\n')
        break
end
% UBC Inversion


end

% figure;plot(misfit)
% figure;plot(2:oo-1,phiWxyz(2:oo-1),2:oo-1,phiWc(2:oo-1),'r')