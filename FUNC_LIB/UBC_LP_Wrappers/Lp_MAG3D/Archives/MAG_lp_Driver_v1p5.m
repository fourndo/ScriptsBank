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
root_lib = 'C:\Users\dominiquef\Dropbox\DOM_Projects\Lp_MAG3D\functions';
addpath (root_lib);

%% USER INPUTS

home_dir = 'C:\Projects\Research\Nigel\Milligan'   ;

% Iteration mode
argin = 'continuous'; 

% Target chi factor
chifact = 1.0;

% Lp-norm on model
q = 0.0:0.1:2.0;

% Lp-norm on gradient (l2-norm by default)
p = 0.0:0.1:2.0;

% Scale phi_cxyy and phi_cm
l = 0.1%0.25:0.25:0.75;

% Trade-off cooling schedule
cool_beta = 0.5;

% Maximum number of iterations
iter_max = 20;

% Treshold values
epsilon = 1e-6;    % Model disturbance
delta = 1e-6;
            
%% DRIVER
% Save driver root

start_message(p,q,l,chifact,cool_beta,iter_max)
 
dos (['mkdir ' home_dir '\Workspace']);

work_dir = [home_dir '\Workspace'];

cd (home_dir);

% model_true = importdata('Model_intrusive.sus');
%% Cycle through all the p, q and lamba 
for ll = 1% : length(l)
    
    for pp = 21% 1 : length(p)

        for qq = 11% 1 : length(q)
            
            
            dos (['copy maginv3d_001.sus' ' ' work_dir]);
            dos (['copy maginv3d.log' ' ' work_dir]);   
            
            
            % Create new directory for current parameters
            [file_list,new_dir] = create_dir(p(pp),q(qq),l(ll),home_dir,root_lib);
            
            % Read all the input files and extract parameters
            [ndata,beta_in,L_scale,phid,phi,meshfile,obs_file,logfile]=read_input;
            
            % Load topocheck file
            load topo_model.txt
%             load dist_weight.txt
            obs = importdata(obs_file, ' ', 3);
            data = obs.data;
            
            [mesh]=get_UBCmesh(meshfile);
            
            dX = mesh(3,mesh(3,:)~=0)';
            dY = mesh(4,mesh(4,:)~=0)';
            dZ = mesh(5,mesh(5,:)~=0)';
            
            nX = mesh(1,1); %size(X,1);    %number of cell in X
            nY = mesh(1,2); %size(X,2);    %number of cell in Y
            nZ = mesh(1,3); %size(X,3);    %number of cell in Z

            cd (work_dir);
            
            
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
            beta = beta_in/2;     % Trade-off parameter

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

            oo=1;
            strnum = '001';
            model = load('maginv3d_001.sus');
            
%% Iteration loop computing weightings
            while oo <= iter_max && phid(oo) >= ndata * chifact;


            % Import the file

                
                model = model .* topo_model;


            % Compute GRAD(m) for all three components
                Wxm = comp_gradxm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wym = comp_gradym(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wzm = comp_gradzm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);

            % Compute dimension matrix for the GRAD terms (FOR DEV ONLY)
                dXx = get_dXx(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dYy = get_dYy(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dZz = get_dZz(model, nX, nY, nZ, dX, dY, dZ ,topo_model);


                phi_xyz = alphax * (dXx.*Wxm)' * ((dXx.*Wxm)) +...
                    alphay * (dYy.*Wym)' * ((dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * ((dZz.*Wzm));
% 
%                 
%                 phi_xyz = alphax * (Wxm)' * ((Wxm)) +...
%                     alphay * (Wym)' * ((Wym)) +...
%                     alphaz * (Wzm)' * ((Wzm));
                
                
                % lp-norm on gradient ||GRAD(m)||p
                wcx = 1 ./ ( (abs(Wxm) ) .^ ( 2 - p(pp) ) + delta );

                
                wcy = 1 ./ ( (abs(Wym) ) .^ ( 2 - p(pp) ) + delta );

                
                wcz = 1 ./ ( (abs(Wzm) ) .^ ( 2 - p(pp) ) + delta );

                phi_cxyz = alphax * (dXx.*Wxm)' * (wcx .* (dXx.*Wxm)) +...
                    alphay * (dYy.*Wym)' * (wcy .* (dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * (wcz .* (dZz.*Wzm));
                
                
%                 phi_cxyz = alphax * (Wxm)' * (wcx .* (Wxm)) +...
%                     alphay * (Wym)' * (wcy .* (Wym)) +...
%                     alphaz * (Wzm)' * (wcz .* (Wzm));
               
                % Scale norm of ||GRAD(m)||p on ||GRAD(m)||2
                scale_xyz = phi_xyz / phi_cxyz;

                
                % Minimum support on model || m ||p
                wcm = 1 ./ ( ( abs(model) ) .^ ( 2 - q(qq) ) + epsilon );

                
                % Compute norms and scales
                phi_cm = Ws'.*model' * (wcm .* (Ws .* Wr.* model));

                scale_cm = (scale_xyz * phi_cxyz) / (phi_cm);
                                
                phi_m = scale_cm*phi_cm + scale_xyz*phi_cxyz;
                phi_iter = phid + beta * phi_m;

                % Final weighting vectors
                wcm = (lambda * scale_cm * wcm);

                wcx = abs(1-lambda) * (scale_xyz * wcx);
                wcy = abs(1-lambda) * (scale_xyz * wcy);
                wcz = abs(1-lambda) * (scale_xyz * wcz);
                
               % Adjust the trade-off parameter for next iteration            
                    if phid(oo) < ndata*2 && oo >1
                      beta =cool_beta *beta;
                    elseif oo >= 1  
                      beta = cool_beta * beta;
                    end


%% Save weighting file

                cd (work_dir);

                save('input_w.dat', 'wcm','-ascii')
                fid=fopen('input_w.dat','a');
                fprintf(fid,'\n');
                fclose(fid);

                save('input_w.dat', 'wcx','-ascii','-append')
                fid=fopen('input_w.dat','a');
                fprintf(fid,'\n');
                fclose(fid);

                save('input_w.dat', 'wcy','-ascii','-append')
                fid=fopen('input_w.dat','a');
                fprintf(fid,'\n');
                fclose(fid);


                save('input_w.dat', 'wcz','-ascii','-append')
                fid=fopen('input_w.dat','a');
                fprintf(fid,'\n');
                fclose(fid);

%% Create input file

                %Update ctrl file
                fid=fopen('inv_ctrl.dat','w');
                fprintf(fid,'2 ! mode\n');
                fprintf(fid,'%12.8f 1.0   ! par tolc\n',beta);
                fprintf(fid,'%s  ! observations file\n',[home_dir '\' obs_file]);
                fprintf(fid,'%s \n',[home_dir '\maginv3d.mtx']);
                fprintf(fid,'%s  ! initial model\n',[work_dir '\maginv3d_' strnum '.sus'] );
                fprintf(fid,'VALUE 0  ! reference model\n');
                fprintf(fid,'NULL  ! active cell file\n');
                fprintf(fid,'VALUE 0  ! lower bounds\n');
                fprintf(fid,'VALUE 1  ! upper bounds\n');
                fprintf(fid,'%d %d %d  ! Le, Ln, Lz\n',Lx,Ly,Lz);
                fprintf(fid,'SMOOTH_MOD\n');
                fprintf(fid,'input_w.dat  ! weighting file');
                fclose(fid);

                oo=oo+1;


%% Call maginv3d_newbounds and invert...

                switch argin 
                    case 'continuous'
                        
                    dos('maginv3d_newbounds inv_ctrl.dat');
                                       
                    [phid(oo),iter]=read_log(logfile);
                    
                    if phid(oo) == 99999
                        
                        sprintf('Inversion has stopped at convergence. Iteration: %i\n',oo)
                        continue
                        
                    elseif phid(oo) == -1
                        
                        phid(oo) = 99999;
                        
                    else
                        
                        
                    % Rename file
                    
                    if iter>=10
                        
                        strnum = ['0' num2str(iter)];
                        
                    else
                        
                        strnum = ['00' num2str(iter)];
                        
                    end
                    
                    dos (['copy maginv3d_' strnum '.sus' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\maginv3d_' strnum '.sus maginv3d_0' num2str(oo) '.sus ']);

                    dos (['copy maginv3d_' strnum '.pre' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\maginv3d_' strnum '.pre maginv3d_0' num2str(oo) '.pre ']);

                    dos (['copy maginv3d.log' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\maginv3d.log maginv3d_0' num2str(oo) '.log ']);
            
% 
                    model = load(['maginv3d_' strnum '.sus']);
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