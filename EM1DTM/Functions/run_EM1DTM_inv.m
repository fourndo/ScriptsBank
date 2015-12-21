function [m_out, m_misfit, phid, phim, beta,pred,bHSpace] = run_EM1DTM_inv(work_dir,meshfile,data,pred,Q,phid,beta,m,m_ref,iter,txloc,txheight,rxloc,rxheight,uncert,component,tc,wf)
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
alpha = [0.1 1]; % alpha_s  alpha_z
tol = []; % convergence test
hankle = []; % Hankle transform parameter
fourier = []; % Fourier transform parameter
output = 2; % amount of output

%% Load 3D mesh
[meshfile]=get_UBC_mesh([work_dir '\' meshfile]);
nx = meshfile(1,1); %size(X,1);    %number of cell in X
ny = meshfile(1,2); %size(X,2);    %number of cell in Y
nz = meshfile(1,3); %size(X,3);    %number of cell in Z

mcell = nx*ny*nz;

% Cell size array
dx = meshfile(3,1:nx);
dy = meshfile(4,1:ny);
dz = meshfile(5,1:nz);

%% Re-shape model
m = reshape(m,nz,nx,ny);
m_ref = reshape(m_ref,nz,nx,ny);

%% Pre-allocate for inversion output
nsnds = size(data,1);

% Create final 1D model result matrix. At most nz cells from 3D mesh
m_out = ones(nz,nx,ny)*1e-8;
m_misfit = ones(nz,nx,ny)*1e-8;
bHSpace     = ones(nz,nx,ny)*1e-8;

% Pre-allocate results
itern   = zeros(nsnds,1);
phim    = zeros(nsnds,1);
phi     = zeros(nsnds,1);

%% write waveform to file
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
% varargin = varargin{1};
% nvarargin = length(varargin);
% switch nvarargin
%     case 0
%         
%         inimodel    = ones(nz,nx,ny)*1e-3;
%         refsmodel   = ones(nz,nx,ny)*1e-3;
%         reffmodel   = ones(nz,nx,ny)*1e-3;
%         
%     case 1
%         
%         inimodel    = reshape(varargin{1}(:),nz,nx,ny);
%         refsmodel   = inimodel;
%         reffmodel   = inimodel;
% 
%     case 2
%         
%         inimodel    = reshape(varargin{1}(:),nz,nx,ny);
%         refsmodel   = reshape(varargin{2}(:),nz,nx,ny);
%         reffmodel   = refsmodel;
%         
%     case 3
%         
%         inimodel    = reshape(varargin{1}(:),nz,nx,ny);
%         refsmodel   = reshape(varargin{2}(:),nz,nx,ny);
%         reffmodel   = reshape(varargin{3}(:),nz,nx,ny);
%         
% end



%% Run all the stations in a loop
for ii = 1 : nsnds
    
    ntcstr = num2str( sum( data(ii,:)<0 ) );
%     ntcstr = num2str(length(tc));
    nrx = size(rxloc,1);
    ntc = length(tc);
    compstr = cell(nrx,1);
    
    if sum( data(ii,:)<0 )==0
        
        fprintf('STOP! no data at %i\n',ii);
        
    end
    
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
                uncert(p,:) = uncert(p,:) * 1e9; 
            case 'by'
                compstr{p} = [' y ',ntcstr,' 4']; % nT for B field data (who did this? must use SI)
                data(p,:) = data(p,:) * 1e9; 
                uncert(p,:) = uncert(p,:) * 1e9; 
            case 'bz'
                compstr{p} = [' z ',ntcstr,' 4']; % nT for B field data (who did this? must use SI)
                data(p,:) = data(p,:) * 1e9; 
                uncert(p,:) = uncert(p,:) * 1e9; 
            case 'dbxdt'
                compstr{p} = [' x ',ntcstr,' 3']; % V for dB/dt field data
            case 'dbydt'
                compstr{p} = [' y ',ntcstr,' 3']; % V for dB/dt field data
            case 'dbzdt'
                compstr{p} = [' z ',ntcstr,' 3']; % V for dB/dt field data
        end
    end
    
    %% write input file
    fid = fopen([work_dir,'\','em1dtm.in'],'wt');
    fprintf(fid,'em1dtm\n');
    fprintf(fid,'em1dtm.obs\n');
    
%     if iter==1
%         fprintf(fid,'inimodel.con\n');
%         fprintf(fid,'refmodel.con\n');        
%     else
        fprintf(fid,'inimodel.con\n');
        fprintf(fid,'refmodel.con\n');
%     end
    fprintf(fid,'NONE\n');
    fprintf(fid,'NONE\n');
    fprintf(fid,'%f  %f  %f  %f  %f\n',huberekblom);
    fprintf(fid,'%f  %f\n',alpha);
    
%     if iter == 1
%         
%         fprintf(fid,'2\n');
%         fprintf(fid,'1 0.5\n');%beta(ii),beta(ii)/2,0.5);
%         fprintf(fid,'50 \n');
%     
%     else
        
        fprintf(fid,'1\n');
        fprintf(fid,'%f %f %f\n',beta(ii)/10,beta(ii)/5,0.5);
        fprintf(fid,'1\n');
        
%     end
        
    if isempty(tol)
        fprintf(fid,'DEFAULT\n');
    else
        fprintf(fid,'%f\n',tol);
    end
    if isempty(hankle)
        fprintf(fid,'DEFAULT\n');
    else
        fprintf(fid,'%d\n',hankle);
    end
    if isempty(fourier)
        fprintf(fid,'DEFAULT\n');
    else
        fprintf(fid,'%d\n',fourier);
    end
    fprintf(fid,'%d\n',output);
    fclose(fid);

    %% write obs file
    fid = fopen([work_dir,'\','em1dtm.obs'],'wt');
    fprintf(fid,'1\n');
    fprintf(fid,'%f %f %f\n',rxloc(1),rxloc(2),txheight(ii));
    fprintf(fid,'%d  ', length(txloc) );

    for jj = 1 : size(txloc,1)
    fprintf(fid,' %f  %f ',txloc(jj,1:2));
    end

    fprintf(fid,'%f\n',-txheight(ii)); % horizontal tx loop, so all vertex co-planar
    fprintf(fid,'em1dtm.wf\n');
    fprintf(fid,'%d  3\n',nrx); % use second for time
%     for p = 1:nrx % loop over rx
    nnztc = data(ii,:)<0; %nntzc(1:3) = 0;
    temp = data(ii,nnztc==1);
    temptc = tc(nnztc==1);
    tempsd = uncert(ii,nnztc==1);
    for p = 1 : nrx
        fprintf(fid,'1  %f  %f  %f %s\n',rxloc(p,1:2),-rxheight(ii),compstr{p});
        for q = 1 : sum(nnztc)
            fprintf(fid,'  %e  1  %e  v  %e\n',temptc(q),temp(q), tempsd(q));% 0.05 * abs(temp(q)) + tempsd(q));
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
    condin = m(Q(ii,3):end,Q(ii,1),Q(ii,2));
    condref = m_ref(Q(ii,3):end,Q(ii,1),Q(ii,2));
    
    ncells = length(dz_layer);
    
    % Check if phid has reached misfit, if no then continue
    if phid(ii) <= sum(nnztc);
        
        m_out(:,Q(ii,1),Q(ii,2)) = m(:,Q(ii,1),Q(ii,2));
        m_misfit(:,Q(ii,1),Q(ii,2)) =  phid(ii);
        
         beta(ii) = beta(ii);
        fprintf('##\nStation %i has al ready reached the target misfit\n##\n',ii)
        continue
    end
    
    fid1 = fopen([work_dir '\inimodel.con'],'wt');
    fprintf(fid1,'%i\n',ncells+1);
    
    fid2 = fopen([work_dir '\refmodel.con'],'wt');
    fprintf(fid2,'%i\n',ncells+1);
    
    
    for jj = 1 : ncells
        
        % If first iteration then right out best fitting half-space
%         if iter == 1
%             fprintf(fid1,'%12.4f\t%12.8e\n',dz_layer(jj),condin(jj));
%             
%         else
            
            fprintf(fid1,'%12.4f\t%12.8e\n',dz_layer(jj),condin(jj));
            
%         end
        
        fprintf(fid2,'%12.4f\t%12.8e\n',dz_layer(jj),condref(jj));
        
        
    end
    
    fprintf(fid1,'%12.4f\t%12.8e\n',0.0,condin(jj));
    fclose(fid1);
    
    fprintf(fid2,'%12.4f\t%12.8e\n',0.0,condref(jj));
    fclose(fid2);

    %% run code
    cd(work_dir);
    system('em1dtm');
    cd(root_dir);
    fid = fopen([work_dir,'\em1dtm.con'],'r');
    tline = fgetl(fid);
    nlayer = str2num(tline);
    
    model = zeros(nlayer,1);

    for p = 1:nlayer
        tline = fgetl(fid);
        temp = sscanf(tline,'%f');
        model(p) = temp(2);
    end
    fclose(fid);

    % Project back to 3D mesh
    m_out(Q(ii,3):end,Q(ii,1),Q(ii,2)) = model(1:end-1);
    
    % Copy top cell all the way up the mesh to avoid 
    % interpolating air cells later
    m_out(1:Q(ii,3),Q(ii,1),Q(ii,2)) = model(1);
    
    % read pred
    fid = fopen([work_dir,'\em1dtm.prd'],'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
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

    % read em1dtm.out
    fid = fopen([work_dir,'\em1dtm.out'],'r');
    line = fgets(fid);
    while line~=-1
        
        iteration = regexp(line,'Iteration','match');
        
        if isempty(iteration) == 0
            
            invout      = regexp(line,'(?<=\=)[a-zA-Z_0-9+.- ][^,]*','match');
            phid(ii)    = str2num(invout{1});
            beta(ii)    = str2num(invout{2});
            phim(ii)    = str2num(invout{3});
%             phi(ii)     = str2num(invout{4});
            
        end
        
        temp = regexp(line,'Best-fitting','match');
        if isempty(temp) == 0
            
            invout      = regexp(line,'\d*\.?\d*','match');
            bHSpace(:,Q(ii,1),Q(ii,2)) = str2double(invout{1})*10^(-str2double(invout{2}));
            
        end
        
        line = fgets(fid);
        
    end
    fclose(fid);

    m_misfit(:,Q(ii,1),Q(ii,2)) =  phid(ii);
%     cd(oldFolder); % get beck to previous directory
end
    
m_out = reshape(m_out,nz*nx*ny,1);

m_misfit = reshape(m_misfit,nz*nx*ny,1);
bHSpace  = reshape(bHSpace,nz*nx*ny,1);
end