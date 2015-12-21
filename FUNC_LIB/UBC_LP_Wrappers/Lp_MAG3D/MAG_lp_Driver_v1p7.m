% MAG_lp_Driver_v1.6
% Author: D Fournier
% Last Update: July 11th, 2013
%
% This program invert for "compact" model using a iterative lp-norm on the
% model and gradient. It is an adapted algorithm based on the classic 
% paper of Last & Kubic. The central solver is the UBC - MAGINV3D_v5 code. 
% Apply compactness by creating a w.dat file, after computing Ws, Wx, Wz
% values. 
%
% //INPUT// Control File
% # home_dir: Directory for the location of the mesh, topo, maginv3d.mtx
%
% # UBC_dir: Directory for the location of the unconstrained models
% 
% # p [VALUE or RANGE between 0:2 ]: 
%   p-norm applied on the gradient: ||GRAD(m)||p (p = 2 recommended for MAG)
% 
% # q [VALUE or RANGE between 0:2 ]: 
%   q-norm applied on the model:||m||q (q = 0 recommended for MAG)
%
% # l [VALUE or RANGE between 0:2 ]: 
%   Ratio between ||GRAD(m)||p and ||m||q. (l = 1 recommended for MAG)
%                           l=1 is equivalant to equal weight on both terms
%
% # cooling schedule between iterations (usually 0.5 for UBC)
% 
% # target chi-factor  (usually 1.0 for UBC)
% 
% # Iter_start: Iteration to start with from the unconstrained inversion 
%   Recommended to start with 1 or 2 for MAG
%
% //OUTPUT//
% (for "step" mode)
% w.dat: Weighting matrix computed
% input.inp: Input file for the maginv3d
%
% (for "continuous" mode)
% maginv3d.sus: Inverted mag models

close all
clear all

root_lib = 'C:\Users\dominiquef\Dropbox\DOM_Projects\Lp_MAG3D\MAG_functions';
addpath (root_lib);

%% USER INPUTS

input_file = 'C:\Projects\3796_AGIC_Research\MAG3D\Scott\lp_MAG3D_input.inp'   ;

%% Read input file
% Iteration mode

[home_dir,UBC_dir,iter_start,chifact,pvec,qvec,lvec,cool_beta,iter_max] = read_input(input_file);

            
%% Initialization
dos (['mkdir ' home_dir '\Workspace']);

work_dir = [home_dir '\Workspace'];

cd (home_dir);

argin = 'continuous'; 

% model_true = importdata('Model_intrusive.sus');
%% Cycle through all the p, q and lamba 
for ll = 1 : length(lvec)
    
    for pp = 1 : length(pvec)

        for qq = 1 : length(qvec)
            
            start_message(pvec(pp),qvec(qq),lvec(ll),chifact,cool_beta,iter_max)
            
            [file_list,new_dir] = create_dir(pvec(pp),qvec(qq),lvec(ll),home_dir,root_lib);
            
            dos (['copy ' UBC_dir '\maginv3d_00' num2str(iter_start) '.sus' ' ' work_dir]);
            dos (['rename ' work_dir '\maginv3d_00' num2str(iter_start) '.sus maginv3d_001.sus ']);
            dos (['copy ' UBC_dir '\maginv3d.log' ' ' work_dir]);   
                 
            % Read all the input files and extract parameters
            [meshfile,obsfile,topofile]=get_input_files;
            
            % Load topocheck file
            load topo_model.txt
            
            % Load mesh and create dimension vectors
            [mesh]=get_UBCmesh(meshfile);
            
            dX = mesh(3,mesh(3,:)~=0)';
            dY = mesh(4,mesh(4,:)~=0)';
            dZ = mesh(5,mesh(5,:)~=0)';
            
            nX = mesh(1,1); %size(X,1);    %number of cell in X
            nY = mesh(1,2); %size(X,2);    %number of cell in Y
            nZ = mesh(1,3); %size(X,3);    %number of cell in Z

            cd (work_dir);
            
            [ndata,beta_in,L_scale,phid,mref,mref_file,iter]=read_log_MAGv4(home_dir,'maginv3d.log',iter_start);
            
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


            lambda = lvec(ll);     % Scaling between smooth and compact
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
            Wr = ones(mcell,1);%dist_weight;

            oo=1;
            strnum = '001';
            model = load('maginv3d_001.sus');
            
            
            % Treshold values
            epsilon = 1e-8;    % Model disturbance
            delta = 1e-8;
            
%% Iteration loop computing weightings
            while oo <= iter_max && phid(oo) >= ndata * chifact;


            % Import the file

                m = model-mref;
                model = model .* topo_model;
                m = m .* topo_model;

            % Compute GRAD(m) for all three components
                Wxm = comp_gradxm(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                Wym = comp_gradym(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                Wzm = comp_gradzm(model, nX, nY, nZ, dX, dY, dZ ,topo_model);

            % Compute dimension matrix for the GRAD terms (FOR DEV ONLY)
                dXx = get_dXx(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dYy = get_dYy(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dZz = get_dZz(model, nX, nY, nZ, dX, dY, dZ ,topo_model);


                phi_xyz = alphax * (dXx.*Wxm)' * ((dXx.*Wxm)) +...
                    alphay * (dYy.*Wym)' * ((dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * ((dZz.*Wzm));
                
                
                % lp-norm on gradient ||GRAD(m)||p
                wcx = 1 ./ ( (abs(Wxm) ) .^ ( 2 - pvec(pp) ) + delta );

                
                wcy = 1 ./ ( (abs(Wym) ) .^ ( 2 - pvec(pp) ) + delta );

                
                wcz = 1 ./ ( (abs(Wzm) ) .^ ( 2 - pvec(pp) ) + delta );

                phi_cxyz = alphax * (dXx.*Wxm)' * (wcx .* (dXx.*Wxm)) +...
                    alphay * (dYy.*Wym)' * (wcy .* (dYy.*Wym)) +...
                    alphaz * (dZz.*Wzm)' * (wcz .* (dZz.*Wzm));
               
               
                % Scale norm of ||GRAD(m)||p on ||GRAD(m)||2
                scale_xyz = phi_xyz / phi_cxyz;

                
                % Minimum support on model || m ||p
                wcm = 1 ./ ( ( abs(m) ) .^ ( 2 - qvec(qq) ) + epsilon );

                
                % Compute norms and scales
                phi_cm = alphas * (Ws'.*m') * (wcm .* (Ws .* m));

                scale_cm = (scale_xyz * phi_cxyz) / (phi_cm);
                                
                phi_m = alphas*scale_cm*phi_cm + scale_xyz*phi_cxyz;
                phi_iter = phid + beta * phi_m;

                % Final weighting vectors
                wcm = (lambda * scale_cm * wcm);

                wcx = abs(2-lambda) * (scale_xyz * wcx);
                wcy = abs(2-lambda) * (scale_xyz * wcy);
                wcz = abs(2-lambda) * (scale_xyz * wcz);
                
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
                fprintf(fid,'%s  ! observations file\n',[home_dir '\' obsfile]);
                fprintf(fid,'%s \n',[home_dir '\maginv3d.mtx']);
                fprintf(fid,'%s  ! initial model\n',[work_dir '\maginv3d_' strnum '.sus'] );
                if length(mref)==1
                    fprintf(fid,'VALUE %8.3f  ! reference model\n',mref);
                else
                    fprintf(fid,'%s  ! reference model\n',mref_file);
                end
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
                                       
                    [arg,arg,arg,phid(oo),arg,arg,iter]=read_log(home_dir,'maginv3d.log','end');
                    clear arg
                    
                    if phid(oo) == 99999
                        
                        sprintf('Inversion has stopped at convergence. Iteration: %i\n',oo)
                       continue
                        
                        
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