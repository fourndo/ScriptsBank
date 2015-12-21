function [ndata,beta_in,L_scale,phid,mref,mref_file,iter]=read_log(home_dir,logfile,argin)
% Opens UBC log file and extract misfit and iteration number.

L_scale=0;
beta_in=0;
ndata=0;
phid = [];
iter = 1;
mref_file = [];
fid=fopen(logfile,'rt');

line=0;
count = 1;
% Go through the log file and extract data and the last achieved misfit
while line~=-1         	
line=fgets(fid); %gets next line 

    if line==-1

        break

    end

    if isempty(findstr(line,'# of data:'))==0
        
            ndata = str2num(line(31:end));
        
    end

    if isempty(findstr(line,'beta:'))==0

            beta_in = str2num(line(21:end));
        
    end

    % Extract alpha values
    if isempty(findstr(line,'Le, Ln, Lz:'))==0

            L_scale = str2num(line(25:end));
        
    end


    if isempty(findstr(line,'   data misfit:'))==0

            phid = str2num(line(21:end));

    end

    if isempty(findstr(line,'Iteration:'))==0
            iter = str2double(line(21:end));
            phid = [];
        
    end
    
    if isempty(findstr(line,'Reference Model:'))==0
                    
            if isempty(findstr(line,'VALUE'))==0
                mref = str2double(line(28:end));
                
            elseif isempty(findstr(line,':\'))==0
                mref_file = strtrim(line(23:end));
                mref = load(mref_file);
            else
                mref_file = strcat(home_dir, strtrim(line(23:end)));
                mref = load(mref_file);
            end
        
    end
    

    
    switch argin 
        case 'end'
            
            continue 
            
        otherwise
                
                if iter==str2num(argin) && isempty(phid)==0
                    break
                end
    end
 count= count+1;
end

    if isempty(phid)==1
        
        phid=99999;
%         fprintf('File ended at line %i\n',count);
%         fprintf('Did not find the information needed - review log file\n')
        
    end

fclose(fid);
