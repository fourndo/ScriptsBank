% MAG_lp_Driver_v1.2
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
%% Library call
close all
clear all
root_lib =  'C:\Users\dominiquef\Dropbox\DOM_Projects\Lp_ZTEM3D\functions';
addpath (root_lib);

%% USER INPUTS

home_dir = 'C:\Projects\3796_AGIC_Research\ZTEM3D'   ;
UBC_dir = 'C:\Projects\3796_AGIC_Research\ZTEM3D\UBC_results' ;

% Iteration mode
argin = 'step';   

% Target chi factor
chifact = 1.0;

% Lp-norm on model
q = 0.0;

% Lp-norm on gradient (l2-norm by default)
p = 2;

% Scale phi_cxyy and phi_cm
l = 0.5;

% Trade-off cooling schedule
cool_beta = 0.2;

% Maximum number of iterations
iter_max = 15;


%% Initialization

start_message(p,q,l,chifact,cool_beta,iter_max)
 
dos (['mkdir ' home_dir '\Workspace']);

work_dir = [home_dir '\Workspace'];

cd (home_dir);

%% Cycle through all the p, q and lamba 
for ll = 1 : length(l)
    
    for pp = 1 : length(p)

        for qq = 1: length(q)

             % Create new directory for current parameters
            [file_list,new_dir] = create_dir(p(pp),q(qq),l(ll),home_dir,root_lib); 
                 
            % Read all the input files and extract parameters
            [meshfile,obsfile,topofile]=read_input;
            
            % Load topocheck file
            load topo_model.txt
            
            [mesh]=get_UBCmesh(meshfile);
            
            dX = mesh(3,mesh(3,:)~=0)';
            dY = mesh(4,mesh(4,:)~=0)';
            dZ = mesh(5,mesh(5,:)~=0)';
            
            nX = mesh(1,1); %size(X,1);    %number of cell in X
            nY = mesh(1,2); %size(X,2);    %number of cell in Y
            nZ = mesh(1,3); %size(X,3);    %number of cell in Z

            cd (work_dir);

                switch argin 
                    case 'continuous'
                        oo=1;
                        iter=1;
                        
                        dos (['del ' home_dir '\Workspace /Q']);
                        dos (['copy ' UBC_dir '\inv_0' num2str(iter_start) '.con' ' ' work_dir]);
                        dos (['rename ' work_dir '\inv_0' num2str(iter_start) '.con' ' inv_01.con ']);
                        dos (['copy ' UBC_dir '\mt3dinv.log' ' ' work_dir]);  
            
                    case 'step'
                        
                        oo=1;
                        iter = input('Inversion in mode "step". Enter current iteration:\n');

                            if iter ==1
                                
                                % Starting model from UBC unconstrained inversion
                                
                                iter_start = input('First iteration. Specify starting iteration from the UBC inversion:\n');
                                
                                dos (['del ' home_dir '\Workspace /Q']);
                                dos (['copy ' UBC_dir '\inv_0' num2str(iter_start) '.con' ' ' work_dir]);
                                dos (['rename ' work_dir '\inv_0' num2str(iter_start) '.con' ' inv_01.con ']);
                                dos (['copy ' UBC_dir '\mt3dinv.log' ' ' work_dir]);  
                                
                            else
                                
                                dos (['del ' home_dir '\Workspace /Q']);
                                dos (['copy ' home_dir '\inv_01.con' ' ' work_dir]);
                                dos (['copy ' home_dir '\mt3dinv.log' ' ' work_dir]);
                                iter_start = 1;
                                
                            end
                                
                    otherwise
                        fprintf('Mode must be "step" or "continuous"\n');
                        break
        
                end

            [ndata,beta_in,alpha,phid,bkgd_model,ref_model]=read_log('mt3dinv.log',iter_start);
            
%% Set-up parameters

            mcell=nX*nY*nZ;      %Number of model cells  

%             % %Length scales
%             Lx=L_scale(1);         
%             Ly=L_scale(2);
%             Lz=L_scale(3);

            % Scaling factors
            alphas = alpha(1);    %Coefficient applied to WstWs in UBC
            alphac = alphas;    % Weight on model minimum support
            alphax = alpha(2);    % Weight on gradx 
            alphay = alpha(3);    % Weight on grady
            alphaz = alpha(4);    % Weight on gradz

            epsilon = 1e-11;    % Model disturbance 
            lambda = l(ll);     % Scaling between smooth and compact
            beta = beta_in;     % Trade-off parameter

            % Create derivative matrices
            mnull = []; % No null values for now

            wcx=ones((nX-1)*nY*nZ,1);
            wcy=ones(nX*(nY-1)*nZ,1);
            wcz=ones(nX*nY*(nZ-1),1);  

            Ws = getWs_3D(mcell,dX,dY,dZ,topo_model);

            Wr = ones(mcell,1);%dist_weight;
            
            strnum = (['0' num2str(iter_start)]);
            model = load('inv_01.con');
%% Iteration loop computing weightings
            while oo <= iter_max %&& phid(oo) >= ndata;

                
            % Import the file

%                 model = load(model_file);
%                 model = model .* topo_model;

%                 model = model - ref_model;
                model = model.*topo_model;
                
            % Compute GRAD(m) for all three components
                Wxm = comp_gradxm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wym = comp_gradym(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);
                Wzm = comp_gradzm(model, nX, nY, nZ, dX, dY, dZ ,topo_model,Wr);

                model = model - ref_model;
                model = model.*topo_model;
                
            % Compute dimension matrix for the GRAD terms
                dXx = get_dXx(model, nX, nY, nZ, dX, dY, dZ , Wr);
                dYy = get_dYy(model, nX, nY, nZ, dX, dY, dZ , Wr);
                dZz = get_dZz(model, nX, nY, nZ, dX, dY, dZ , Wr);

                phi_xyz = alphax * (dXx .* Wxm)' * Wxm +...
                    alphay * (dYy .* Wym)' * Wym +...
                    alphaz * (dZz .* Wzm)' * Wzm;

                % Minimum support on gradient

                wcx = 1 ./ ( abs(Wxm) .^ ( 2 - p(pp) ) + epsilon );

                wcy = 1 ./ ( abs(Wym) .^ ( 2 - p(pp) ) + epsilon );

                wcz = 1 ./ ( abs(Wzm) .^ ( 2 - p(pp) ) + epsilon );

                phi_cxyz = alphax * (dXx .* Wxm)' * (wcx .* Wxm) +...
                    alphay * (dYy .* Wym)' * (wcy .* Wym) +...
                    alphaz * (dZz .* Wzm)' * (wcz .* Wzm);

                scale_xyz = phi_xyz / phi_cxyz;

                % Minimum support on model
                wcm = 1 ./ ( abs( model ) .^ ( 2 - q(qq) ) + epsilon );

                % Compute norms and scales
                phi_cm = alphas * (Ws .* model )' * (wcm .* (Ws .* model) );

                scale_cm = (scale_xyz * phi_cxyz) / (phi_cm);

                phi_m = scale_cm*alphac*phi_cm + scale_xyz*phi_cxyz;

                phi_iter = phid + beta * phi_m;

                % Final weighting vectors
                wcm = (lambda * scale_cm * wcm);

                wcx = abs(1-lambda) * ( scale_xyz * wcx);
                wcy = abs(1-lambda) * ( scale_xyz * wcy);
                wcz = abs(1-lambda) * ( scale_xyz * wcz);

               % Adjust the trade-off parameter            
                    if phid(oo) < ndata*2 && oo*iter >1
                      beta = 0.75*beta;
                    elseif oo*iter>1  
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
                fid=fopen('mt3dinv.inp','w');
                fprintf(fid,'%s ! mesh file\n',(meshfile));
                fprintf(fid,'%10.8f ! background conductivity file\n',bkgd_model);
                fprintf(fid,'%s ! observations file\n',(obsfile));
                fprintf(fid,'inv_01.con ! initial model file\n');
                fprintf(fid,'%10.8f ! reference model file\n',ref_model);
                fprintf(fid,'TOPO_FILE %s ! topography (active cell file for adjusted topo)\n',topofile);
                fprintf(fid,'BOUNDS_NONE ! BOUNDS_NONE | BOUNDS_CONST bl bu | BOUNDS_FILE file\n');
                fprintf(fid,'input_w.dat  ! weighting file\n');           
                fprintf(fid,'%6.3f %6.3f %6.3f ! tradeoff parameter |(max) (min) (cooling multiplier)|\n',beta,beta,cool_beta);
                fprintf(fid,'%6.3f %6.3f %6.3f %6.3f ! alpha_s  alpha_x alpha_y alpha_z\n',alphas,alphax,alphay,alphaz);
                fprintf(fid,'NOT_CHANGE_MREF\n');
                fprintf(fid,'SMOOTH_MOD\n');
                fprintf(fid,'USE_LOG_COND  ! model type\n');
                fprintf(fid,'1  ! chifact\n');
                fprintf(fid,'0  ! method\n');
                fprintf(fid,'0.01  0.001  3  ! tol_nl, mindm, nit\n');
                fprintf(fid,'0.01  3  ! intol, max_linit\n');
                fprintf(fid,'5e-008  5e-008  ! fortol, inintol\n');
                fprintf(fid,'5  0.01  0.01  ! max_it_bicg, droptol, droptolWTW\n');
                fclose(fid);

                oo=oo+1;


%% Call maginv3d_newbounds and invert...

                switch argin 
                    case 'continuous'
                        
                    dos('maginv3d_newbounds inv_ctrl.dat');
                    [ndata,beta_in,L_scale,phid(oo)]=read_log('dcinv3d.log',1);
                    dos (['copy inv_01.con' ' ' home_dir '\' new_dir]);
                    dos (['rename ' home_dir '\' new_dir '\inv_01.con inv_00' num2str(oo) '.con ']);

%                     dos (['copy maginv3d_001.pre' ' ' home_dir '\' new_dir]);
%                     dos (['rename ' home_dir '\' new_dir '\maginv3d_001.pre maginv3d_00' num2str(oo) '.pre ']);
% 
%                     dos (['copy maginv3d.log' ' ' home_dir '\' new_dir]);
%                     dos (['rename ' home_dir '\' new_dir '\maginv3d.log maginv3d_00' num2str(oo) '.log ']);

                    model = load(['inv_' strnum '.con']);

%                     residual(pp,qq,ll) = norm((model - model_true.*topo_model),'inf');
%                     finalmodels(pp,qq,ll,:) = model;
                
%% Extract misfit after inversion
                    [ndata,beta_in,alpha,phid,back_model,ref_model,bkgd_model,ref_model]=read_log('mt3dinv.log',1);

                   case 'step'
                       
                    break
                    
                end
                
                fclose('all');

                if oo>iter_max
                    fprintf('Inversion has reached the maximum number of iterations: %i\n',iter_max)
                    fprintf('Verify convergence\n')
                    fprintf('**End of inversion for lp: %i , lq: %i , lambda %4.2e**\n', p(pp), q(qq) , lambda)
                end


            end

        end

    end
    
end

% save('MAG3D_model_Intrusive.mat','finalmodels','l','p','q','phid');