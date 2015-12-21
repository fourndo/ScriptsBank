function [ndata,beta_in,alphas,phid,mref_file,iter,topofile]=read_log_DC2D_v5(home_dir,work_dir,logfile,argin)
% Opens dcinv2d.log file and extract parameters.

cd(work_dir)

alphas=[];
beta_in=0;
ndata=0;
phid = [];
iter = 1;
mref_file = [];
fid=fopen(logfile,'rt');
topofile =[];

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

    if isempty(findstr(line,'chosen beta ='))==0

            beta_in = str2num(line(20:end));
        
    end

    % Extract alpha values
    if isempty(findstr(line,'ALPHA VALUE'))==0

            alphas = str2num(line(12:end));
        
    end


    if isempty(findstr(line,'  achieved misfit ='))==0

            phid = str2num(line(20:end));

    end

    if isempty(findstr(line,'Iteration '))==0
            iter = str2double(line(11:end));
            phid = [];
        
    end
    
    if isempty(findstr(line,'REF_MOD'))==0
                    
            if isempty(findstr(line,'VALUE'))==0
                mref_file = str2double(line(14:end));
                
            elseif isempty(findstr(line,'DEFAULT'))==0
                mref_file = 0;
                
            elseif isempty(findstr(line,'FILE'))==0
                mref_file = strcat(home_dir, strtrim(line(13:end)));
                
                
            end
        
    end
    
    if isempty(findstr(line,'TOPO'))==0
          
        if isempty(findstr(line,'DEFAULT'))==0
            
            topofile = 'DEFAULT';
            
        elseif isempty(findstr(line,'FILE'))==0
            
           topofile = ['FILE ' home_dir '\' strtrim(line(11:end))];
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
