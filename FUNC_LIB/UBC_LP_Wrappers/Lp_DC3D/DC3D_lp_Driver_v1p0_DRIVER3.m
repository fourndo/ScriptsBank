% MAG_lp_Driver_v1.2
% Author: D Fournier
% Last Update: July 2th, 2013
%
% This program invert for "compact" model using a iterative lp-norm on the
% model and gradient. It is an adapted algorithm based on the classic 
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

close all
clear all
root_lib = 'C:\Users\dominiquef\Dropbox\DOM_Projects\Lp_DC3D\functions';
addpath (root_lib);

%% USER INPUTS

home_dir = 'C:\Projects\3796_AGIC_Research\DCIP3D'   ;
UBC_dir = 'C:\Projects\3796_AGIC_Research\DCIP3D\UBC_Results' ;

% Iteration mode
argin = 'continuous'; 

% Iteration from the unconstrained inversion to use
iter_start = 1;

% Depth weighting
wdepth = [16;8;4;2;1];


ini_model = 1e-3;
ref_model = 1e-3;

% Target chi factor
chifact = 0.01;

% Lp-norm on gradient (l2-norm by default)
p = 2;

% Lq-norm on model
q = 0;

% Scale phi_cxyy and phi_cm
l = 1.5;%0.25:0.25:0.75;

% Trade-off cooling schedule
cool_beta = 0.5;

% Maximum number of iterations
iter_max = 20;

% Treshold values
epsilon = 1e-8;    % Model disturbance
delta = 1e-8;
     


%% DRIVER
% Save driver root

start_message(p,q,l,chifact,cool_beta,iter_max)
 
dos (['mkdir ' home_dir '\Workspace3']);

work_dir = [home_dir '\Workspace3'];

cd (home_dir);

% model_true = importdata('Model_intrusive.sus');
%% Cycle through all the p, q and lamba 
for ll = 1 : length(l)
    
    for pp = 1 : length(p)

        for qq = 1 : length(q)
            
            % Create new directory for current parameters
            [file_list,new_dir] = create_dir(p(pp),q(qq),l(ll),home_dir,root_lib);
            
            dos (['del ' home_dir '\Workspace /Q']);
            dos (['copy ' UBC_dir '\dcinv3d_0' num2str(iter_start) '.con' ' ' work_dir]);
%             dos (['rename ' work_dir '\dcinv3d_0' num2str(iter_start) '.con dcinv3d_01.con ']);
            dos (['copy ' UBC_dir '\dcinv3d.log' ' ' work_dir]);   
                 
            % Read all the input files and extract parameters
            [meshfile,obsfile,topofile]=read_input;
            
            % Load topocheck file
            load topo_model.txt
%             load dist_weight.txt
%             obs = importdata(obs_file, ' ', 3);
%             data = obs.data;
            
            [mesh]=get_UBCmesh(meshfile);
            
            dX = mesh(3,mesh(3,:)~=0)';
            dY = mesh(4,mesh(4,:)~=0)';
            dZ = mesh(5,mesh(5,:)~=0)';
            
            nX = mesh(1,1); %size(X,1);    %number of cell in X
            nY = mesh(1,2); %size(X,2);    %number of cell in Y
            nZ = mesh(1,3); %size(X,3);    %number of cell in Z

            cd (work_dir);
            
            [ndata,beta_in,L_scale,phid]=read_log('dcinv3d.log',iter_start);
            
%% Set-up parameters

            mcell=nX*nY*nZ;      %Number of model cells  

            % %Length scales
            Lx=L_scale(1);         
            Ly=L_scale(2);
            Lz=L_scale(3);

            % Scaling factors
            alphas = 1/Lx^2;    %Coefficient applied to WstWs in UBC
            alphac = alphas;    % Weight on model minimum support
            alphax = 1.0e-0;    % Weight on gradx 
            alphay = 1.0e-0;    % Weight on grady
            alphaz = 1.0e-0;    % Weight on gradz


            lambda = l(ll);     % Scaling between smooth and compact
            beta = beta_in;     % Trade-off parameter

            % Create derivative matrices
            mnull = []; % No null values for now

            wcx=ones((nX-1)*nY*nZ,1);
            wcy=ones(nX*(nY-1)*nZ,1);
            wcz=ones(nX*nY*(nZ-1),1);  

%             Wx = getWx_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
%             Wy = getWy_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
%             Wz = getWz_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
            
            Ws = getWs_3D(mcell,dX,dY,dZ,topo_model);
            
            % Get depth weighting function
%             Wr = get_Wr(mesh,data,topo_model);
            Wr = ones(mcell,1);%dist_weight;

            % Create depth weighting matrices
            
            [wxdepth,wydepth] = get_depthweight(wdepth,nZ,nX,nY,topo_model);
            
            oo=1;
            strnum = (['0' num2str(iter_start)]);
            model = load(['dcinv3d_' strnum '.con']);
            
%% Iteration loop computing weightings
            while oo <= iter_max && phid(oo) >= ndata * chifact;


            % Import the file

%                 model = 1./(abs(log10(model)));
%                 model = model - ref_model;


            % Compute GRAD(m) for all three components
                Wxm = comp_gradxm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wym = comp_gradym(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wzm = comp_gradzm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);

            % Compute dimension matrix for the GRAD terms (FOR DEV ONLY)
                dXx = get_dXx(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dYy = get_dYy(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dZz = get_dZz(model, nX, nY, nZ, dX, dY, dZ ,topo_model);


                phi_xyz = alphax * (wxdepth .* dXx.*Wxm)' * ((dXx.*Wxm)) +...
                    alphay * (wydepth .* dYy.*Wym)' * ((dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * ((dZz.*Wzm));
                
                
                % lp-norm on gradient ||GRAD(m)||p
                wcx = wxdepth ./ ( (abs(Wxm) ) .^ ( 2 - p(pp) ) + delta ) ;

                
                wcy = wydepth ./ ( (abs(Wym) ) .^ ( 2 - p(pp) ) + delta );

                
                wcz = 1 ./ ( (abs(Wzm) ) .^ ( 2 - p(pp) ) + delta );

                phi_cxyz = alphax * (dXx.*Wxm)' * (wcx .* (dXx.*Wxm)) +...
                    alphay * (dYy.*Wym)' * (wcy .* (dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * (wcz .* (dZz.*Wzm));
                
               
                % Scale norm of ||GRAD(m)||p on ||GRAD(m)||2
                scale_xyz = phi_xyz / phi_cxyz;

                
                % Minimum support on model || m ||p
                wcm = 1 ./ ( ( abs(model) ) .^ ( 2 - q(qq) ) + epsilon );

                
                % Compute norms and scales
                phi_cm = Ws'.*model' * (wcm .* (Ws .* Wr.* model));

%                 scale_cm = (Ws'.*model' * ((Ws .* model))) / (phi_cm);
                        
                    scale_cm = (scale_xyz * phi_cxyz) / (alphas *phi_cm);
                    
                phi_m = scale_cm*phi_cm + scale_xyz*phi_cxyz;
                
                phi_iter = phid + beta * phi_m;

                % Final weighting vectors
                wcm = (lambda * scale_cm * wcm);

                wcx = abs(2-lambda) * (scale_xyz * wcx) ;
                wcy = abs(2-lambda) * (scale_xyz * wcy) ;
                wcz = abs(2-lambda) * (scale_xyz * wcz);
%                 wcm = (scale_cm * wcm);
% 
%                 wcx = (scale_xyz * wcx) ;
%                 wcy = (scale_xyz * wcy) ;
%                 wcz = (scale_xyz * wcz);
 
                
                
               % Adjust the trade-off parameter for next iteration            
                    if phid(oo) < ndata*2 && oo >1
                      beta =cool_beta *beta;
                    elseif oo >= 1  
                      beta = cool_beta * beta;
                    end


%% Save weighting file

                cd (work_dir);

                save('input_w.dat', 'wcm','-ascii')
%                 fid=fopen('input_w.dat','a');
%                 fprintf(fid,'\n');
%                 fclose(fid);

                save('input_w.dat', 'wcx','-ascii','-append')
%                 fid=fopen('input_w.dat','a');
%                 fprintf(fid,'\n');
%                 fclose(fid);

                save('input_w.dat', 'wcy','-ascii','-append')
%                 fid=fopen('input_w.dat','a');
%                 fprintf(fid,'\n');
%                 fclose(fid);


                save('input_w.dat', 'wcz','-ascii','-append')
%                 fid=fopen('input_w.dat','a');
%                 fprintf(fid,'\n');
%                 fclose(fid);

%% Create input file

                %Update ctrl file
                fid=fopen('inv_ctrl.dat','w');
                fprintf(fid,'1 0 !!! restart=1, max # of iterations\n');
                fprintf(fid,'2, %12.8f !!! mode chifact \n',beta);
                fprintf(fid,'%s  ! observations file\n',[home_dir '\' obsfile]);
                fprintf(fid,'%s  !!! mesh or ncel,aspr\n',[home_dir '\' meshfile]);
                fprintf(fid,'%s !!! topography (NULL = no topography)\n',[home_dir '\' topofile]);
                fprintf(fid,'%s !!! initial model (NULL = same as reference model)\n',['dcinv3d_' strnum '.con']);
                fprintf(fid,'%10.8f !!! reference model (NULL = calculate)\n',ref_model);
                fprintf(fid,'null !!! active cells file (NULL = all cells are active)\n');
                fprintf(fid,'%d %d %d  ! Le, Ln, Lz\n',Lx,Ly,Lz);
                fprintf(fid,'daub2 !!! sensitivity wavelet type: Daubechies-4, 2 vanishing moments\n');
                fprintf(fid,'1 0.02 !!! relative constrcution error applied\n');
                fprintf(fid,'input_w.dat  ! weighting file\n');
                fprintf(fid,'0 !!! store sensitivity matrix in disk\n');
                fprintf(fid,'1e-005 !!! tolerance\n');
                fprintf(fid,'-1 !!! no limit to use memeory\n');
                fclose(fid);

                oo=oo+1;


%% Call maginv3d_newbounds and invert...

                switch argin 
                    case 'continuous'
                        
                    dos('dcinv3d inv_ctrl.dat');
                                       
                    [ndata,beta_in,L_scale,phid(oo)]=read_log('dcinv3d.log',1);
                    
                    if phid(oo) == 99999
                        
                        sprintf('Inversion has stopped at convergence. Iteration: %i\n',oo)
                        continue
                        
                    elseif phid(oo) == -1
                        
                        phid(oo) = 99999;
                        
                    else
                        
                        
                    % Rename file
                    strnum = '01';
                    
                    dos (['copy dcinv3d_01.con' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\dcinv3d_01.con dcinv3d_0' num2str(oo) '.con ']);

                    dos (['copy dcinv3d.pre' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\dcinv3d.pre dcinv3d_0' num2str(oo) '.pre ']);

                    dos (['copy dcinv3d.log' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\dcinv3d.log dcinv3d_0' num2str(oo) '.log ']);
            
% 
                    model = load(['dcinv3d_' strnum '.con']);
%                     dos (['copy maginv3d_001.sus' ' ' home_dir '\' new_dir]);
%                     dos (['rename ' home_dir '\' new_dir '\maginv3d_001' strnum '.sus maginv3d_0' num2str(oo) '.sus ']);

%                     model = load(['maginv3d_001.sus']);
                                       
                    end
                    
%                     misfit(pp,qq,ll) = phid(end);
%                     finalmodels(pp,qq,ll,:) = model;
                
                    % Extract misfit after inversion

                    

                   case 'step-by-step'

                    otherwise

                        fprintf('Input argument not recognize\n')
                        fprintf('Should be "step-by-step" or "continuous"\n')
                        break
                end
                
%                 if mag3d_iter >= 6
%                     break
%                 end
                
                fclose('all');

                if oo>iter_max
                    fprintf('Inversion has reached the maximum number of iterations: %i\n',iter_max)
                    fprintf('Verify convergence\n')
                    fprintf('**End of inversion for lp: %i , lq: %i , lambda %4.2e**\n', p(pp), q(qq) , lambda)
                end


            end
                    cd ..
            
        end

    end
    
end

% save('MAG3D_model_Intrusive_v2.mat','finalmodels','l','p','q','phid');