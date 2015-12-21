function [meshfile,obsfile,topofile,mstart,mref,chi_target,alphas,beta,nlayer,FLAG] = MAG3C_OMES_read_inp(inputfile)
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
        
        topofile = strtrim(arg{1});
        
        if isempty(regexp(topofile,'null','match'))==0 || isempty(regexp(topofile,'NULL','match'))==0
            topofile = [];
            fprintf('No topography - entire mesh used.\n');
        else
            fprintf('Topography: \t\t %s\n',topofile);
        end
        
        
    elseif count == 4
        
        mstart = strtrim(arg{1});
        
        if isempty(regexp(mstart,'VALUE','match'))==0
            
            temp = regexp(mstart,'\s','split');
            mstart = str2num(temp{2});
            fprintf('Starting model: \t\t %e\n',mstart);
            
        else
            
            fprintf('Starting model: \t\t %s\n',mstart);
            
        end
        
        
    elseif count == 5
        
        mref = strtrim(arg{1});
        
        if isempty(regexp(mref,'VALUE','match'))==0
            
            temp = regexp(mref,'\s','split');
            mref = str2num(temp{2});
            fprintf('Starting model: \t\t %e\n',mref);
            
        else
            
            fprintf('Starting model: \t\t %s\n',mref);
            
        end
        
    
    elseif count == 6
        
        chi_target = str2num(arg{1});
        
        if isempty(chi_target)==1
            
            chi_target = 1.0;
            
        end
        
        fprintf('Target chi factor: \t %f\n',chi_target);
        
    elseif count == 7
        
        alphas = str2num(arg{1});
        
        if isempty(alphas)==1 || length(alphas)~=4
            
            fprintf('Error in input file at line %i\n',count);
            fprintf('Requires four numerical values (e.g.--> 0.001 1 1 1\n')
            break
            
        end
        
        fprintf('Alpha values: \t\t %4.1e %4.1e %4.1e %4.1e\n',alphas);
        
    elseif count == 8
        
       beta = str2num(arg{1});   
       
       if isempty(beta)==1
            
          fprintf('Starting Beta: \t\t COMPUTED\n'); 
          
       else
           
          fprintf('Starting Beta: \t\t %f\n',beta); 
       end
        
        
    elseif count == 9

        nlayer = strtrim(arg{1});
        
        if isempty(regexp(nlayer,'DEFAULT','match'))==0
            
            nalyer = 1;
            fprintf('Number of OMES layers set to default: 1\n');
       
        elseif isempty(str2num(arg{1}))==0 && length(str2num(arg{1}))==1
            
            nlayer = str2num(arg{1});
            
            fprintf('Number of OMES layers: \t\t\t %i',nlayer)

            fprintf('\n');
            
            

        else
            
            nlayer = str2num(arg{1});
            fprintf('Error in input file at line %i\n',count);
            fprintf('Require single numerical value | DEFAULT (1)\n')
            break
            
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