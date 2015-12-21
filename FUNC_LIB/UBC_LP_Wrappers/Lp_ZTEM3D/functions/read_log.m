function [ndata,beta_in,alpha,phid,bkgd_model,ref_model]=read_log(logfile,iter_start)
% Opens UBC log file and extract misfit and iteration number.

alpha = [];
beta_in = [];
ndata = [];
phid = [];
iter = 1;

fid=fopen(logfile,'rt');


max_num_lines = 30000;
% Go through the log file and extract data and the last achieved misfit
for ii=1:max_num_lines         	
line=fgets(fid); %gets next line 

    if line==-1
        fprintf('File ended at line %i\n',ii);
        fprintf('Did not find the information needed - review log file\n')
        break
    end

    
    if length(strtrim(line))>=length('Iteration:')
        description = strtrim(line);
        if strcmp(description(1:10),'Iteration:')==1
            
            iter = str2double(description(11:end));
            phid = [];
        end
    end
    
    if length(strtrim(line))>=length('# of data:')
        description = strtrim(line);
        if strcmp(description(1:10),'# of data:')==1
            ndata = str2num(description(11:end));
        end
    end

    if length(strtrim(line))>=length('========== BETA = ')
        description = strtrim(line);
        if strcmp(description(1:18),'========== BETA = ')==1
            beta_in = str2num(description(19:31));
            
            iter = str2double(description(34:38));
        end
        
    end

    % Extract alpha values
    if length(strtrim(line))>=length('alpha (s,x,y,z):')
        description = strtrim(line);
        if strcmp(description(1:16),'alpha (s,x,y,z):')==1
            alpha = str2num(description(18:end));
        end
    end


    if length(strtrim(line))>=length('phi_d:')
        description = strtrim(line);
        if strcmp(description(1:6),'phi_d:')==1
            phid = str2num(description(7:end));
        end
    end
    
    if length(strtrim(line))>=length('reference model was set to a constant of:')
        description = strtrim(line);
        if strcmp(description(1:41),'reference model was set to a constant of:')==1
            ref_model = str2num(description(42:end));
        end
    end
    
    if length(strtrim(line))>=length('reference model was read from file:')
        description = strtrim(line);
        if strcmp(description(1:35),'reference model was read from file:')==1
            ref_model = strtrim(description(36:end));
        end
    end
    
    if length(strtrim(line))>=length('background model was set to a constant of:')
        description = strtrim(line);
        if strcmp(description(1:42),'background model was set to a constant of:')==1
            bkgd_model = str2num(description(43:end));
        end
    end

    
    if iter==iter_start && isempty(phid)==0
        break
    end

end

fclose(fid);
