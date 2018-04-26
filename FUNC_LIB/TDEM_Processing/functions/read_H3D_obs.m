function [dobs,tx,index] = read_H3D_obs(obsfile)
% Read obs and predicted data and write to file 
% INPUT           
% work_dir: directory where to find data.obs and data.pre
% obsfile: name of observation data file
%
% OUTPUT
% X,Y,Z,Time, ..., dBx/dt, uncert, dBy/dt, uncert , dBz/dt, uncert
% Repeat for dBx, dBy, dBz


rfid=fopen(obsfile,'rt');

line=fgets(rfid);

while isempty(strfind(line,'N_TRX'))==1
          line=fgets(rfid);
end

ntrx = str2num(line(7:end));

tx = zeros(ntrx,6);
count = 0;
for oo = 1 : ntrx

%     while isempty(strfind(line,'TRX_LOOP'))==1
%           line=fgets(rfid);
%           
%           
%     end
%     
%     while isempty(str2num(line))==1
%               line=fgets(rfid);
%               
%     end
%     
%     tx(oo,:) = str2num(line);

    while isempty(strfind(line,'N_RECV'))==1
          line=fgets(rfid);
    end

    nrecv = str2num(line(7:end));
    line=fgets(rfid);
    ntime = str2num(line(7:end));

    for ii = 1 : nrecv

        count = count + 1;
        for jj = 1 : ntime
            
            line=fgets(rfid); 
            
            while length(line)<2
                line=fgets(rfid);
            end        

            if oo==1 && ii==1 && jj==1

                index = count;
                dobs = str2num(line);

            else

                index = [index;count];
                dobs = [dobs;str2num(line)];

            end
                  


        end

    end

end
fclose(rfid);

end

