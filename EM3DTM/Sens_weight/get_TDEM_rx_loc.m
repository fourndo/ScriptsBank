% Function get_TDEM_rx_loc(filename)

clear all 
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\4264_Vedanta_TEM\Inversion\50m';
obsfile = '4264_Composite_Data_Zcomp.dat';

fid = fopen([work_dir '\' obsfile],'r');
line = fgets(fid);

while isempty(regexp(line,'N_TRX','match'))==1;
    
    line = fgets(fid);
    
end
ntrx = str2double(regexp(line,'\d*','match'));

for ii = 1 : ntrx
    
    while isempty(regexp(line,'N_RECV','match'))==1;
    
    line = fgets(fid);
    
    end
    

    nrec = str2double(regexp(line,'\d*','match'));
    line = fgets(fid);
    ntimes = str2double(regexp(line,'\d*','match'));
    line= fgets(fid);
    countrx = 1;
    for jj = 1 : nrec

        for kk = 1 : ntimes
            
            % Skip blank lines
           while isempty(regexp(line,'\d*','match'));
               
               line = fgets(fid);
               
           end
           
           if kk == 1
               
               data = str2num(line);
               rx{ii}(countrx,1:3) = data(1:3);
               countrx = countrx+1;
               
           end
        
           line = fgets(fid);
           
        end
        
    end
    
    
    
    
end 