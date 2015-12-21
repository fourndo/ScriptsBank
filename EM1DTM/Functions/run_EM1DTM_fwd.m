function [pred] = run_EM1DTM_fwd(work_dir,meshfile,data,Q,txloc,txheight,rxloc,rxheight,sd,component,tc,wf,varargin)
% function [model, pred, phid, beta, phim,] = em1dtm(txloc,txheight,rxloc,rxheight,data,sd,component,tc,wf,layers,varargin)
% call tem iD inversion code em1dtm as a function
% this is very low-level; does not support object
% all the variables are numeric, string, etc.
%
%                 [  point1x  point1y  ]
%         txloc = [  point2x  point2y  ]
%                 [    ...      ...    ]
%                 [  pointnx  pointny  ]
%
%      txheight = [ TxHeight ] 
%
%                 [    RX1x     RX1y   ]
%         rxloc = [    RX2x     RX2y   ]
%                 [    ...      ...    ]
%                 [    RXnx     RXny   ]
%
%      rxheight = [ RxHeight ] 
%
%                 [ DataRX1TC1  DataRX1TC2  ...  DataRX1TCn ]
%          data = [ DataRX2TC1  DataRX2TC2  ...  DataRX2TCn ]
%                 [     ...        ...      ...       ...   ]
%                 [ DataRXnTC1  DataRXnTC2  ...  DataRXnTCn ]
%
%                 [ StdDevRX1TC1  StdDevRX1TC2  ...  StdDevRX1TCn ]
%            sd = [ StdDevRX2TC1  StdDevRX2TC2  ...  StdDevRX2TCn ]
%                 [      ...          ...       ...         ...   ]
%                 [ StdDevRXnTC1  StdDevRXnTC2  ...  StdDevRXnTCn ]
%
%     component = { RX1Comp1; RX2Comp2; ... RXnCompn},
%              where RX?Comp? can be 'Bx'/'By'/'Bz'/'dBxdt'/'dBydt'/'dBzdt'
%
%            tc = [TC1, TC2, ..., TCn]
%
%            wf = [Time1, Amp1; Time2, Amp2; ...; Timek, Ampk] or 'STEPOFF'
%
%        layers = [thick1; thick2; ...; thickp]
%                   NOTE: a 0 will be added to the end of layers to
%                   represent the basement in the inversion, but the
%                   conductivity of basement will not be returned.
%
%    OPTIONAL ARGUMENTS:
%
%      inimodel = [con1; con2; ...; conp]
%                   if missing, use best fitting halfspace
%
%     refsmodel = [con1; con2; ...; conp]  % Mref for smallest component
%                   if missing, the same as inimodel
%
%     reffmodel = [con1; con2; ...; conp]  % Mref for flattest component
%                   if missing, the same as refsmodel
%
% 
% ORIGINAL CODE FROM: SEOGI KANG
% ADAPTED BY: Dom Fournier (2014-03-21)
% 
root_dir = pwd; % get the absolute path of this file
% oldFolder = cd(thisfile_path); % get into modeling directory

% internal parameters
huberekblom = [1000  2 0.0001  2 0.0001]; % Huber and Ekblom parameters
alpha = [0.0001 1]; % alpha_s  alpha_z
tol = []; % convergence test
hankle = []; % Hankle transform parameter
fourier = []; % Fourier transform parameter
output = 2; % amount of output

%% Load 3D mesh
[xn,yn,zn]=read_UBC_mesh([work_dir '\' meshfile]);

nx = length(xn) - 1 ;
ny = length(yn) - 1 ;
nz = length(zn) - 1 ;

dz = zn(1:end-1) - zn(2:end);
%% Pre-allocate for inversion output
nsnds = size(data,1);

% Pre-allocate results
pred = zeros(size(data,1),size(data,2)); % same size as data and sd

%% write wf file
fid = fopen([work_dir '\em1dtm.wf'],'wt');
if strcmpi(wf,'STEPOFF')
    fprintf(fid,'STE\n');
else
    fprintf(fid,'%d\n',size(wf,1));
    for p = 1:size(wf,1)
        fprintf(fid,'%15.5e   %15.5e\n',wf(p,:));
    end
end
fclose(fid);

%% Load starting, reference model if provided
varargin = varargin{1};
nvarargin = length(varargin);
switch nvarargin
    case 0
        
        inimodel    = ones(nz,nx,ny)*1e-4;
        refsmodel   = ones(nz,nx,ny)*1e-4;
        reffmodel   = ones(nz,nx,ny)*1e-4;
        
    case 1
        
        inimodel    = reshape(varargin{1}(:),nz,nx,ny);
        refsmodel   = inimodel;
        reffmodel   = inimodel;

    case 2
        
        inimodel    = reshape(varargin{1}(:),nz,nx,ny);
        refsmodel   = reshape(varargin{2}(:),nz,nx,ny);
        reffmodel   = refsmodel;
        
    case 3
        
        inimodel    = reshape(varargin{1}(:),nz,nx,ny);
        refsmodel   = reshape(varargin{2}(:),nz,nx,ny);
        reffmodel   = reshape(varargin{3}(:),nz,nx,ny);
        
end

%% Run all the stations in a loop
for ii = 1 : nsnds
    
    ntcstr = num2str(length(tc));
    nrx = size(rxloc,1);
    ntc = length(tc);
    compstr = cell(nrx,1);
    if ischar(component)
        component = {component};
    end

    compstr = cell(nrx,1);

    if ischar(component)
        component = {component};
    end

    for p = 1:nrx
        switch lower(component{p})
            case 'bx'
                compstr{p} = [' x ',ntcstr,' 4']; % nT for B field data (who did this? must use SI)
                data(p,:) = data(p,:) * 1e9; 
                sd(p,:) = sd(p,:) * 1e9; 
            case 'by'
                compstr{p} = [' y ',ntcstr,' 4']; % nT for B field data (who did this? must use SI)
                data(p,:) = data(p,:) * 1e9; 
                sd(p,:) = sd(p,:) * 1e9; 
            case 'bz'
                compstr{p} = [' z ',ntcstr,' 4']; % nT for B field data (who did this? must use SI)
                data(p,:) = data(p,:) * 1e9; 
                sd(p,:) = sd(p,:) * 1e9; 
            case 'dbxdt'
                compstr{p} = [' x ',ntcstr,' 3']; % V for dB/dt field data
            case 'dbydt'
                compstr{p} = [' y ',ntcstr,' 3']; % V for dB/dt field data
            case 'dbzdt'
                compstr{p} = [' z ',ntcstr,' 3']; % V for dB/dt field data
        end
    end
    
    %% write inp file
    fid = fopen([work_dir,'\','em1dtmfwd.in'],'wt');
    fprintf(fid,'em1dtm.obs\n');
    fprintf(fid,'inimodel.con\n');
    fprintf(fid,'DEFAULT\n');
    fprintf(fid,'DEFAULT\n');
    fprintf(fid,'NO\n');
    fclose(fid);

    %% write obs file
    fid = fopen([work_dir,'\','em1dtm.obs'],'wt');
    fprintf(fid,'%d  ', length(txloc) );

    for jj = 1 : size(txloc,1)
    fprintf(fid,' %f  %f ',txloc(jj,1:2));
    end

    fprintf(fid,'%f\n',-txheight(ii)); % horizontal tx loop, so all vertex co-planar
    fprintf(fid,'em1dtm.wf\n');
    fprintf(fid,'%d  3\n',nrx); % use second for time
%     for p = 1:nrx % loop over rx
    nnztc = ones(size(data,2),1); %nntzc(1:3) = 0;
    temp = data(ii,nnztc==1);
    temptc = tc(nnztc==1);
    tempsd = sd(nnztc==1);
    for p = 1 : nrx
        fprintf(fid,'1  %f  %f  %f %s\n',rxloc(p,1:2),-rxheight(ii),compstr{p});
        for q = 1 : sum(nnztc)
            fprintf(fid,'%e  1 \n',temptc(q));
        end
    end
%     for p = 1 : nrx
%         fprintf(fid,'1  %f  %f  %f z %f 4\n',rxloc(p,1:2),-rxheight(ii),13);
%         for q = 3 : size(data,2)
%             fprintf(fid,'  %e  1  %e  v  %e\n',tc(q),data(ii,q),sd(q));
%         end
%     end
    fclose(fid);

    %%
    % Create layer conductivity, susc and weights
    
    dz_layer = dz(Q(ii,3):end)';
    condin = inimodel(Q(ii,3):end,Q(ii,1),Q(ii,2));
    condref = refsmodel(Q(ii,3):end,Q(ii,1),Q(ii,2));
    
    ncells = length(dz_layer);
    
    fid1 = fopen([work_dir '\inimodel.con'],'wt');
    fprintf(fid1,'%i\n',ncells+1);
    
    
    
    for jj = 1 : ncells
        
        % If first iteration then right out best fitting half-space
       
        fprintf(fid1,'%12.4f\t%12.8e\n',dz_layer(jj),condin(jj));
        
        
    end
    
    fprintf(fid1,'%12.4f\t%12.8e\n',0.0,condin(jj));
    fclose(fid1);
    

    %% run code
    cd(work_dir);
    system('em1dtmfwd');
    cd(root_dir);
    
    % read pred
    fid = fopen([work_dir,'\em1dtmfwd.out'],'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    tline = fgetl(fid);
    
    for p = 1:nrx
        rxline = fgetl(fid);
        temp = regexp(rxline,'\s*','split');
        npred = str2num(temp{end-1});
        
        for q = 1:npred
            tline = fgetl(fid);
            temp = sscanf(tline,'%f');
            pred(ii,q) = temp(3);
        end
        
        % Move time channels to proper slot
        pred(ii,nnztc==1) = pred(ii,1:npred);
        pred(ii,nnztc==0) = nan;
%         if strcmp(rxline(end),'4') % nT for B field
%             pred(ii,1:npred) = pred(ii,1:npred); % go back to SI
%         end
    end
    fclose(fid);

%     cd(oldFolder); % get beck to previous directory
end
