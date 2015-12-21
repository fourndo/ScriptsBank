function [ndata,beta_in,L_scale,phid]=read_log(logfile,iter_start)
% Opens UBC log file and extract misfit and iteration number.

L_scale = 0;
beta_in = 0;
ndata   = 0;
phid    = [];
iter    = 1;

fid=fopen(logfile,'rt');


max_num_lines = 30000;
% Go through the log file and extract data and the last achieved misfit
for ii=1:max_num_lines         	
line=fgets(fid); %gets next line 

    if line==-1
        fprintf('File ended at line %i\n',ii);
        fprintf('Did not find the information needed - review log file\n')
        phid = 99999;
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

    if length(strtrim(line))>=length('tradeoff par:')
        description = strtrim(line);
        if strcmp(description(1:13),'tradeoff par:')==1
            beta_in = str2num(description(14:end));
        end
    end

    % Extract alpha values
    if length(strtrim(line))>=length('Le, Ln, Lz:')
        description = strtrim(line);
        if strcmp(description(1:11),'Le, Ln, Lz:')==1
            L_scale = str2num(description(12:end));
        end
    end


    if length(strtrim(line))>=length('achieved misfit:')
        description = strtrim(line);
        if strcmp(description(1:16),'achieved misfit:')==1
            phid = str2num(description(17:end));
        end
    end

    
    if iter==iter_start && isempty(phid)==0
        break
    end

end

fclose(fid);
