function write_EM1DTM_obs(work_dir,filename,data,tc,txloc,xyz,uncert)
% WRITEEM1DFMOBS(rootname,data,varargin)
% Creates observation file in UBC format for the 1D inversion of frequency
% domain airborn data. 
% 
% INPUTS:
%
% data: all time channels
% xyz: location of observation
% uncert: assigned uncertainties
% loop: 
% 
% 
%% Write obs for all stations\
nstn = size(data,1);

fid = fopen([work_dir '\' filename],'wt');
fprintf(fid,'%i\n',nstn);
   
std_floor = 0.01*max(abs(data),[],1);

for ii = 1 : nstn

    
    fprintf(fid,'%f %f %f\n',xyz(ii,1),xyz(ii,2),-xyz(ii,4));
    
    fprintf(fid,'%d  ', size(txloc,1) );

    for jj = 1 : size(txloc,1)
    fprintf(fid,' %f  %f ',txloc(jj,1) + xyz(ii,1), txloc(jj,2) + xyz(ii,2));
    end

    fprintf(fid,'%f\n',-xyz(ii,4)); % horizontal tx loop, so all vertex co-planar
    fprintf(fid,'em1dtm.wf\n');
    fprintf(fid,'1  3\n'); % use second for time
%     for p = 1:nrx % loop over rx
    nnztc = data(ii,:)<0; %nntzc(1:3) = 0;
    temp = data(ii,nnztc==1);
    temptc = tc(nnztc==1);
    tempsd = uncert(ii,nnztc==1) + std_floor(nnztc==1);

    fprintf(fid,'1  %f  %f  %f z %i 3\n',xyz(ii,1),xyz(ii,2),-xyz(ii,4),sum(nnztc));
    for q = 1 : sum(nnztc)
        fprintf(fid,'  %e  1  %e  v  %e\n',temptc(q),temp(q), tempsd(q));% 0.05 * abs(temp(q)) + tempsd(q));
    end

%     fprintf(fid,'\n');
    
end

fclose(fid);