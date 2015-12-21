%% Format Downhole DCIP data and write to UBC format
% Written by: Dom Fourner
% Last update: 2013-12-08

close all
clear all

data_dir = 'C:\Projects\4160_Abitibi_Windfall\Data\';
work_dir = 'C:\Projects\4160_Abitibi_Windfall\Processing\';
datafile = '14-cols-input.dat';

% Load data line by line
fid = fopen([data_dir datafile],'rt');
line=fgets(fid);
line=fgets(fid);

% Pre-allocate space for data array: Atlas
% Transmitter location <X1,Y1,Z1,X2,Y2,Z2>
Atlas.tx = zeros(1,6);

countid = 0;
counttx = 1;
while (line)~=-1
    
        
    % Check for beggining of block
    if line(1)=='/'
        
        holeid = regexp(line(3:end),' - ','split');
        
        P1 = char(strtrim(holeid(1)));
        P2 = char(strtrim(holeid(2)));
                
        countA = 1;     % Count along hole data in A
        countB = 1;     % Count along hole data in B
        countAB = 1;    % Count across hole data between A-B
        countid = countid+1;
        
    else
        
        data = str2num(line) ;
        
         % Check for transmitter location, stored if first time encountered.
        if sum(Atlas.tx(:,1)==data(1))==0
            
            % Store Tx1(XYZ) , Tx2(XYZ)
            Atlas.tx(counttx,:) = data(1:6);
            txid = counttx;
            
            counttx = counttx+1;
                        
        else
            
            % If already exist, then get the Tx id#
            txid = find(Atlas.tx(:,1)==data(1));

        end
        
                
        if countA==1;
            % Get approximate coordinates of first hole
            DH_A = data(7:9);
        end
        
        % Data is sent to different branchs of Atlas.data depending on:
        % Branch 1: Along hole A -> tx{1} -> DataID -> HoleID -> Data
        %                                 
        %                        -> tx{2} -> DataID -> HoleID -> Data
        %
        % Branch 2: Along hole B ...
        %
        % Branch 3: Cross hole AB ...
        %
        % Compute horizontal distance between poles
        r = ( (DH_A(1) - data(10))^2 + (DH_A(2) - data(11))^2 )^0.5;
        
        % Check if second pole is on the other hole
        % If yes, then could be downhole measurements.
        if r > 50
            
            % Compute distance from both poles.
            r = ( (data(7) - data(10))^2 + (data(8) - data(11))^2 )^0.5;
            
            % If distance is greater than 50, it is an across hole data
            if r > 50
                
                Atlas.data{3}{txid}{countid,1} = ['From_' P1 '_to_' P2];
                Atlas.data{3}{txid}{countid,2}(countAB,:) = data;
                countAB = countAB+1;
                
            else
                
                Atlas.data{2}{txid}{countid,1} = ['Along_' P2 '_with_' P1];
                Atlas.data{2}{txid}{countid,2}(countB,:) = data;
                countB = countB+1;
                
            end
            
            % Else it is a downhole along Hole A
        else
            
            Atlas.data{1}{txid}{countid,1} = ['Along_' P1 '_with_' P2];
            Atlas.data{1}{txid}{countid,2}(countA,:) = data;
            countA = countA+1;
            
        end
        
        
        
        
    end
    
    line=fgets(fid);
           
end

fclose(fid);

%% Write to file in UBC 3D format

% Run through Atlas and write data to file in UBC format.
% Create directory for:
% Along hole P1: Tx1 and Tx2
% Along hole P2: Tx1 and Tx2
% Across hole: Tx1 and Tx2

out_dir = 'C:\Projects\4160_Abitibi_Windfall\Processing';
cross_dir = [out_dir '\AcrossHole'];
down_dir = [out_dir '\DownHole'];
tx_dir{1} = '\Tx1';
tx_dir{2} = '\Tx2';


% Loop over along and across data
for ii = 1 : size(Atlas.data,2)
    
    % Loop over tx's
    for jj = 1 : size(Atlas.data{ii},2)
            
            
        for kk = 1 : size(Atlas.data{ii}{jj},1)

            if isempty(Atlas.data{ii}{jj}{kk,1})==1
                
                continue
                
            end
            
            % If along holes
            if ii==1 || ii ==2
                
                wr2_dir = [down_dir tx_dir{jj}];
                
            else
                
                wr2_dir = [cross_dir tx_dir{jj}];
                
            end
            
            wrt2file = fopen([wr2_dir '\' Atlas.data{ii}{jj}{kk,1} '.dat'],'w');
            fprintf(wrt2file,'!!MIRA - Crosshole DCIP3D BATCH\n');
            %         fprintf(wrt2file,[dc_obs_file '\n']);
            fprintf(wrt2file,'\n\n');
            
            ndata = size(Atlas.data{ii}{jj}{kk,2},1);
            % Write transmiter coordinates
            fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %i\n',...
                Atlas.data{ii}{jj}{kk,2}(1,1),...
                Atlas.data{ii}{jj}{kk,2}(1,2),...
                Atlas.data{ii}{jj}{kk,2}(1,3),...
                Atlas.data{ii}{jj}{kk,2}(1,4),...
                Atlas.data{ii}{jj}{kk,2}(1,6),...
                Atlas.data{ii}{jj}{kk,2}(1,6),...
                ndata);


            
            for ww = 1 : ndata
                
                fprintf(wrt2file,'%15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e %15.8e\n',...
                Atlas.data{ii}{jj}{kk,2}(ww,7),...
                Atlas.data{ii}{jj}{kk,2}(ww,8),...
                Atlas.data{ii}{jj}{kk,2}(ww,9),...
                Atlas.data{ii}{jj}{kk,2}(ww,10),...
                Atlas.data{ii}{jj}{kk,2}(ww,11),...
                Atlas.data{ii}{jj}{kk,2}(ww,12),...
                Atlas.data{ii}{jj}{kk,2}(ww,15),...
                Atlas.data{ii}{jj}{kk,2}(ww,16));
                
                
            end

            fclose(wrt2file);


        end
        
        
    end
    
end

        
        

%Assign coordinates to DC data
%line,station -> line,station,X,Y,Z
