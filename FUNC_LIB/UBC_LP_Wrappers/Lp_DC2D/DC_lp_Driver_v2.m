% DC_lp_Driver_v1.0
% Author: D Fournier
% Last Update: July 2th, 2013
%
% This program invert for "compact" model using a iterative lp-norm on the
% model and gradient. It is an adapted algorithm based on the classic 
% paper of Last & Kubic. The central solver is the UBC - DCINV2D code. 
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
root_lib = 'C:\Users\dominiquef\Dropbox\DOM_Projects\Lp_DC2D\functions';
addpath (root_lib);

%% USER INPUTS

home_dir = 'C:\Users\dominiquef\Dropbox\EOSC556\Chapter2\DCIP2D_Inversion_Exercise\FieldExample'   ;
UBC_dir = 'C:\Users\dominiquef\Dropbox\EOSC556\Chapter2\DCIP2D_Inversion_Exercise\FieldExample\UBC_results';

iter_start='04';

% Target chi factor
chifact = 1.0;

% Lp-norm on model
q = 0:0.1:2;

% Lp-norm on gradient (l2-norm by default)
p = 0:0.1:2;

% Scale phi_cxyy and phi_cm
l = 0.75;%0.0:0.5:1;

% Trade-off cooling schedule
cool_beta = 0.5;

% Maximum number of iterations
iter_max = 40;

% Treshold values
epsilon = 1e-8;    % Model disturbance
delta = 1e-8;

% No-data value for DC
ndval = 1e-8;
            
%% DRIVER
% Save driver root

start_message(p,q,l,chifact,cool_beta,iter_max)
 
dos (['mkdir ' home_dir '\Workspace']);

work_dir = [home_dir '\Workspace'];

cd (home_dir);

% model_true = importdata('Model_intrusive.sus');
%% Cycle through all the p, q and lamba 
for ll = 1% : length(l)
    
    for pp = 21% : length(p)

        for qq = 11%: length(q)

            dos (['copy ' UBC_dir '\dcinv2d_' iter_start '.con' ' ' work_dir '\dcinv2d.con']);
%             dos (['rename ' work_dir 'dcinv2d_' iter_start '.con dcinv2d_01.con']);
            dos (['copy ' UBC_dir '\dcinv2d.log' ' ' work_dir]);   
            dos (['copy ' UBC_dir '\dcinv2d.out' ' ' work_dir]);
            
            
            
            % Create new directory for current parameters
            [file_list,new_dir] = create_dir(p(pp),q(qq),l(ll),home_dir,root_lib);
            
            % Find mesh and obs file in home_directory
            [meshfile,obsfile]=fetch_input(home_dir);
            

            [ndata,beta_in,alpha,phid,mref_file,iter,topofile]=read_log_DC2D_v5(home_dir,work_dir,'dcinv2d.log',iter_start);

            cd(home_dir)
            [mesh]=get_UBCmesh2D(meshfile);
            
            dX = mesh(3,mesh(3,:)~=0)';
%             dY = mesh(4,mesh(4,:)~=0)';
            dZ = mesh(4,mesh(4,:)~=0)';
            
            nX = length(dX); %size(X,1);    %number of cell in X
%             nY = mesh(1,2); %size(X,2);    %number of cell in Y
            nZ = length(dZ); %size(X,3);    %number of cell in Z
            mcell=nX*nZ;      %Number of model cells 
            
            model = loadUBC2D([work_dir '\dcinv2d.con']);
            
            
            % Create topocheck file from input model
            topo = ones(mcell,1);
            topo(model==ndval) = 0;
            
            mref = ones(mcell,1) * mref_file;
            
            
%% Set-up parameters

             cd (work_dir);

            % Scaling factors
            alphas = alpha(1);    %Coefficient applied to WstWs in UBC
            alphac = alphas;    % Weight on model minimum support
            alphax = alpha(2);    % Weight on gradx 
            alphaz = alpha(3);    % Weight on gradz


            lambda = l(ll);     % Scaling between smooth and compact
%             beta = beta_in;     % Trade-off parameter

            % Create derivative matrices
            mnull = []; % No null values for now

            wcx=ones((nX-1)*nZ,1);
%             wcy=ones(nX*(nY-1)*nZ,1);
            wcz=ones(nX*(nZ-1),1);  

%             Wx = getWx_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
%             Wy = getWy_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
%             Wz = getWz_3D(mcell,dX,dY,dZ,reshape(topo_model,nZ,nX,nY));
            
            Ws = getWs_2D(mcell,dX,dZ,topo);
            
            % Get depth weighting function
%             Wr = get_Wr(mesh,data,topo_model);
            Wr = ones(mcell,1);%dist_weight;

            oo=1;

%% Iteration loop computing weightings
            while oo <= iter_max && phid(oo) >= ndata * chifact;


            % Import the file

                
                model(topo==1) = log(model(topo==1)) - log(mref(topo==1)) ;
                model = model.*topo;
%                 model(topo==1) = log(model(topo==1));

            % Compute GRAD(m) for all three components
                Wxm = comp_gradxm2D(model, nX, nZ, dX, dZ ,topo,Wr);
%                 Wym = comp_gradym_v2(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wzm = comp_gradzm2D(model, nX, nZ, dX, dZ ,topo,Wr);

            % Compute dimension matrix for the GRAD terms
                dXx = get_dXx2D(model, nX, nZ, dX, dZ ,topo);
%                 dYy = get_dYy(model, nX, nY, nZ, dX, dY, dZ ,topo_model);
                dZz = get_dZz2D(model, nX, nZ, dX, dZ ,topo);


                phi_xyz = alphax * (dXx.* Wxm)' * (dXx.* Wxm) +...                   
                    alphaz * (dZz .* Wzm)' * (dZz .* Wzm);

                % lp-norm on gradient ||GRAD(m)||p
                wcx = 1 ./ ( (abs(Wxm) ) .^ ( 2 - p(pp) ) + delta );
              
                wcz = 1 ./ ( (abs(Wzm) ) .^ ( 2 - p(pp) ) + delta );

                
                phi_cxyz = alphax * (dXx.* Wxm)' * (wcx .* (dXx.* Wxm)) +...
                           alphaz * (dZz .* Wzm)' * (wcz .* (dZz .* Wzm));
               
                % Scale norm of ||GRAD(m)||p on ||GRAD(m)||2
                scale_xyz = phi_xyz / phi_cxyz;

                
                % Minimum support on model || m ||p
                wcm = 1 ./ ( ( abs(model) ) .^ ( 2 - q(qq) ) + epsilon );

                
                % Compute norms and scales
                phi_cm = Ws'.*model' * (wcm .* (Ws .* Wr.* model));

                scale_cm = (scale_xyz * phi_cxyz) / (phi_cm);
%                 scale_cm = (Ws .* model )' * ((Ws .* model) ) / (phi_cm);
                                
                phi_m = scale_cm*phi_cm + scale_xyz*phi_cxyz;
%                 phi_iter = phid + beta * phi_m;

                % Final weighting vectors
                wcm = (lambda * scale_cm * wcm);
                wcm = reshape(wcm,nZ,nX);
                
                wcx = abs(1-lambda) * (scale_xyz * wcx);
                wcx = reshape(wcx,nZ,nX-1);
%                 wcx = [wcx ones(nZ,1)*mean(mean(wcx))];
                
                wcz = abs(1-lambda) * (scale_xyz * wcz);
                wcz = reshape(wcz,nZ-1,nX);
%                 wcz = [wcz;ones(1,nX)*mean(mean(wcz))];
                
               % Adjust the trade-off parameter for next iteration            
%                     if phid(oo) < ndata*2 && oo >1
%                       beta = 0.75*beta;
%                     elseif oo >= 1  
%                       beta = cool_beta * beta;
%                     end


%% Save weighting file

                cd (work_dir);

                fid=fopen('input_w.dat','w');
                fprintf(fid,'%i %i\n',nX,nZ);
                fclose(fid);
                
                save('input_w.dat', 'wcm','-ascii','-append')
                fid=fopen('input_w.dat','a');
                fprintf(fid,'\n');
                fclose(fid);

                save('input_w.dat', 'wcx','-ascii','-append')
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
                fprintf(fid,'NITER 1\n');
                fprintf(fid,'MESH FILE %s\n',[home_dir '\' meshfile]);
                fprintf(fid,'OBS LOC_X %s\n',[home_dir '\' obsfile]);
                fprintf(fid,'CHIFACT 1\n');
                fprintf(fid,'TOPO %s\n',topofile);
                fprintf(fid,'INIT_MOD FILE dcinv2d.con\n');
                fprintf(fid,'REF_MOD VALUE 1e-3\n');
                fprintf(fid,'ALPHA VALUE 1e-5 1 1\n');
                fprintf(fid,'WEIGHT input_w.dat\n');
                fprintf(fid,'WAVE 2.5e-4 1.0 13\n');
                fprintf(fid,'STORE_ALL_MODELS TRUE\n');
                fprintf(fid,'INVMODE CG\n');
                fprintf(fid,'!CG_PARAM 10 1e-2\n');
                fprintf(fid,'!ACTIVE_CELLS activity.txt\n');
                fprintf(fid,'!HUBER 1e+100\n');
                fprintf(fid,'!EKBLOM 2. 2. 2. 0 0 0\n');
                fprintf(fid,'USE_MREF TRUE\n');
                fprintf(fid,'!BOUNDS VALUE 0.0001 10.0\n');
                fclose(fid);
%                 if isempty(topofile)==1
% 
%                     fprintf(fid,'null  ! topography\n');               
% 
%                 else
% 
%                     fprintf(fid,'%s  ! topography\n',[home_dir '\'  topofile] );
% 
%                 end

                

                oo=oo+1;


%% Call maginv3d_newbounds and invert...
                     
                        
                    dos ('C:\Users\dominiquef\Dropbox\EOSC556\Chapter2\DCIP2D_Inversion_Exercise\dcinv2d inv_ctrl.dat');
                    
                   % Extract misfit after inversion
                    [ndata,beta_in,alpha,phid(oo),mref_file,iter]=read_log_DC2D_v5(home_dir,work_dir,'dcinv2d.log','end');
                       
                       
                    if iter>=10
                        
                        strnum = ['0' num2str(iter)];
                        
                    else
                        
                        strnum = ['00' num2str(iter)];
                        
                    end
                    
                    dos (['copy dcinv2d.con' ' ' home_dir '\' new_dir '\dcinv2d_0' num2str(oo) '.con ']);
%                     dos (['rename ' home_dir '\' new_dir ]);

                    dos (['copy dcinv2d.pre' ' ' home_dir '\' new_dir '\dcinv2d_0' num2str(oo) '.pre ']);
%                     dos (['rename ' home_dir '\' new_dir ]);

                    dos (['copy dcinv2d.log' ' ' home_dir '\' new_dir '\dcinv2d_0' num2str(oo) '.log ']);
%                     dos (['rename ' home_dir '\' new_dir ]);

                                                          
                    
                    
                    model = loadUBC2D('dcinv2d.con');
                    
%                     figure(1) ;imagesc(reshape(model,nZ,nX));
%                     figure(2) ;imagesc(wcm);
%                     residual(pp,qq,ll) = norm((model - model_true.*topo_model),'inf');
%                     finalmodels(pp,qq,ll,:) = model;
                
         
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