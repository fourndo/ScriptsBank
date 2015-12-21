function rx = load_EM1DFM_obs(work_dir,obsfile)
% Load observation file for UBC EM1DFM


%% Load Horizontal dipole experiment

file = [work_dir '\' obsfile];
fid=fopen(file,'rt');    

line=fgets(fid); %gets next line
nsnd = str2num(line); % number of soundings

% Pre-allocate memory
Xsnd  = ones(nsnd,1);
Ysnd  = ones(nsnd,1);
nfreq = ones(nsnd,1);

% Loop over all the soundings
for ss = 1:nsnd
    % Read next line
    temp = str2num(fgets(fid));
    Xsnd(ss) = temp(1);    % X coordinate of transmiters
    Ysnd(ss) = temp(2);    % X coordinate of transmiters
    nfreq(ss) = temp(3);   % Number of frequencies per soundings
    
    if ss==1
        
        freq = zeros(1,nfreq(ss));

        rx = zeros(nsnd,nfreq(ss),1,1,6); %Receiver array

    end
    
    % Loop over all frequencies
    for ff = 1:nfreq(ss)

        line=fgets(fid); %gets next line
        temp = str2num(line);

        freq(ff) = temp(1);
        ntx = temp(2);

        % Loop over all transmitters
        % Pre-allocate memory
        tx_dpm = ones(ntx,1);       %Dipole moment of transmitter
        Z_tx = ones(ntx,1);         %Z location (negative up)

        for tt = 1:ntx

            line=fgets(fid); %gets next line
            temp = regexp(line,'\s+','split');
            temp = temp(1:end);

            tx_dpm(tt) = str2num(temp{1});  
            Z_tx(tt) = str2num(temp{2});     
            nrx = str2num(temp{4});

            % Loop over all receivers
            % Pre-allocate memory
                        
            for rr = 1:nrx
                line=fgets(fid); %gets next line
                temp = regexp(line,'\s+','split');
                temp = temp(1:end);
                rx(ss,ff,tt,rr,1) = str2num(temp{1});   %Dipole momentof rx
                rx(ss,ff,tt,rr,2) = str2num(temp{2});   % x-seperation from tx
                rx(ss,ff,tt,rr,3) = str2num(temp{3});   % y-seperation from tx
                rx(ss,ff,tt,rr,4) = str2num(temp{4});   % z of rx
                rx(ss,ff,tt,rr,5) = str2num(temp{8});   % inphase
                rx(ss,ff,tt,rr,6) = str2num(temp{9});   % quadrature


            end

        end
    end
end
