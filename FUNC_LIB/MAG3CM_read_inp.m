function [meshfile,obsfile,suscfile,chi_target,alphas,beta,pvec,qvec,lvec,FLAG] = MAG3CM_read_inp(inputfile)
% Function [work_dir,meshfile,obsfile,wr_flag,chi_target,alphas,beta,pvec,qvec,lvec] = MAG3D_read_inp('inputfile')
% Read input file for inversion

fid = fopen(inputfile,'r');
line = fgets(fid);

fprintf('Start reading input file\n')
count = 1;
while line~=-1
    
    arg = regexp(line,'!','split');
    
    if count == 1
        
        meshfile = strtrim(arg{1});
        fprintf('Mesh file: \t\t\t %s\n',meshfile);
        
    elseif count == 2
        
        obsfile = strtrim(arg{1});
        fprintf('Observations: \t\t %s\n',obsfile);
      
    elseif count == 3

    suscfile = strtrim(arg{1});
    fprintf('Effective susc: \t\t %s\n',suscfile);
        
        
    elseif count == 4
        
        chi_target = str2num(arg{1});
        
        if isempty(chi_target)==1
            
            chi_target = 1.0;
            
        end
        
        fprintf('Target chi factor: \t %f\n',chi_target);
        
    elseif count == 5
        
        alphas = str2num(arg{1});
        
        if isempty(alphas)==1 || length(alphas)~=4
            
            fprintf('Error in input file at line %i\n',count);
            fprintf('Requires four numerical values (e.g.--> 0.001 1 1 1\n')
            break
            
        end
        
        fprintf('Alpha values: \t\t %4.1e %4.1e %4.1e %4.1e\n',alphas);
        
    elseif count == 6
        
       beta = str2num(arg{1});   
       
       if isempty(beta)==1
            
          fprintf('Starting Beta: \t\t COMPUTED\n'); 
          
       else
           
          fprintf('Starting Beta: \t\t %f\n',beta); 
       end
        
    elseif count == 7
        
        pvec = str2num(arg{1});
        
        if isempty(pvec)==1
            
            fprintf('Error in input file at line %i\n',count);
            fprintf('Require at least one numerical value (e.g.--> 2 )\n')
            break
            
        else
            
            fprintf('lp-norms: \t\t\t ')
            for ii = 1 : length(pvec)
                
                fprintf('%4.2f ',pvec(ii));
                
            end
            fprintf('\n');
            
        end
        
    elseif count == 8

        qvec = str2num(arg{1});

        if isempty(qvec)==1

            fprintf('Error in input file at line %i\n',count);
            fprintf('Require at least one numerical value (e.g.--> 2 )\n')
            break

        else
            
            fprintf('lq-norms: \t\t\t ')
            for ii = 1 : length(qvec)
                
                fprintf('%4.2f ',qvec(ii));
                
            end
            fprintf('\n');
            
        end
        
    elseif count == 9

        lvec = str2num(arg{1});

        if isempty(lvec)==1

            fprintf('Error in input file at line %i\n',count);
            fprintf('Require at least one numerical value (e.g.--> 2 )\n')
            break

        else
            
            fprintf('scaling: \t\t\t ')
            for ii = 1 : length(lvec)
                
                fprintf('%4.2f ',lvec(ii));
                
            end
            fprintf('\n');
            
        end
        
    elseif count == 10

        FLAG = strtrim(arg{1});
        
        if strcmp(FLAG,'SMOOTH_MOD')==0 && strcmp(FLAG,'SMOOTH_MOD_DIF')==0
            
            fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
            fprintf('Please revise input\n')
            return
        
        end
        fprintf('Gradient: \t\t\t %s\n',FLAG);
        
    end
    
    line = fgets(fid);
    count = count +1;
    
end
       
fclose(fid);