function [meshfile,obsfile,topofile,mstart,mref,magfile,weightfile,chi_target,alphas,beta,bounds,norm_vec,eps_p,eps_q,FLAG1,FLAG2,ROT] = MAG3Cinv_read_inp(inputfile)
% Function [work_dir,meshfile,obsfile,wr_flag,chi_target,alphas,beta,pvec,qvec,lvec] = MAG3D_read_inp('inputfile')
% Read input file for inversion

fid = fopen(inputfile,'r');
line = fgets(fid);

fprintf('Start reading input file\n')
count = 1;
while count<17
    
    arg = regexp(line,'!','split');
    
    if count == 1
        
        meshfile = strtrim(arg{1});
        fprintf('Mesh file: \t\t\t %s\n',meshfile);
        
    elseif count == 2
        
        obsfile = strtrim(arg{1});
        fprintf('Observations: \t\t %s\n',obsfile);
    
    elseif count == 3
        
        topofile = strtrim(arg{1});
        
        if isempty(regexp(topofile,'null','match'))==0
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
            fprintf('Reference model: \t\t %e\n',mref);
            
        else
            
            fprintf('Reference model: \t\t %s\n',mref);
            
        end
        
    elseif count == 6
        
         magfile = strtrim(arg{1});
        if isempty(regexp(magfile,'DEFAULT','match'))==0
            magfile = [];
            fprintf('No magnetization file - assume all induced.\n');
        else
            fprintf('Magnetization model: \t\t %s\n',magfile);
        end
        
    elseif count == 7
        
         weightfile = strtrim(arg{1});
        if isempty(regexp(weightfile,'DEFAULT','match'))==0
            weightfile = [];
            fprintf('Cell-based weights - DEFAULT.\n');
        else
            fprintf('Cell-based weights: \t\t %s\n',weightfile);
        end
        
    elseif count == 8
        
        chi_target = str2num(arg{1});
        
        if isempty(chi_target)==1
            
            chi_target = 1.0;
            
        end
        
        fprintf('Target chi factor: \t %f\n',chi_target);
        
    elseif count == 9
        
        alphas = str2num(arg{1});
        
        if isempty(alphas)==1 || length(alphas)~=4
            
            fprintf('Error in input file at line %i\n',count);
            fprintf('Requires four numerical values (e.g.--> 0.001 1 1 1\n')
            break
            
        end
        
        fprintf('Alpha values: \t\t %4.1e %4.1e %4.1e %4.1e\n',alphas);
        
    elseif count == 10
        
       beta = str2num(arg{1});   
       
       if isempty(beta)==1
            
          fprintf('Starting Beta: \t\t COMPUTED\n'); 
          
       else
           
          fprintf('Starting Beta: \t\t %f\n',beta); 
       end
    
   elseif count == 11
        
        temp = strtrim(arg{1});
        
        if isempty(regexp(temp ,'VALUE','match'))==0
            
            temp = regexp(temp,'\s','split');
            bounds(1) = str2num(temp{2});
            bounds(2) = str2num(temp{3});
            fprintf('Bounds: [ %f %f] \n',bounds(1),bounds(2));
            
        else
            
            fprintf('Default bounds: [-1 -1]\n');
            bounds =[-1 1];
            
        end
        
   
    elseif count == 12
        
        if isempty(regexp(arg{1},'VALUE','match'))==0
            
       	    temp = regexp(arg{1},'\s','split');
            for ii = 2 : 6


                    norm_vec(ii-1) = str2num(temp{ii});

                    if isempty(norm_vec(ii-1))==1

                        fprintf('Error in input file at line %i\n',count);
                        fprintf('Requires five numerical value (e.g.--> p ,qx, qy, qz , r )\n')
                        break

                    end



            end
            
            fprintf('Lp-norm: %2.1f %2.1f %2.1f %2.1f %2.1f\n',norm_vec)
        
        elseif isempty(regexp(arg{1},'FILE','match'))==0
            
            temp = regexp(arg{1},'\s','split');
            norm_vec = strtrim(temp{2});
            
            fprintf('Lp-norm from file: %s\n',norm_vec)

            	
        end
        
    elseif count == 13
               
       if isempty(regexp(arg{1},'DEFAULT','match'))==0
          eps_p = [];
          eps_q = [];
          fprintf('Threshold p and q: \t\t COMPUTED\n'); 
          
       elseif isempty(regexp(arg{1},'VALUE','match'))==0
            temp = regexp(arg{1},'\s','split');
            eps_p = str2num(temp{2});   
            eps_q = str2num(temp{3});
            fprintf('Threshold p and q: \t\t %e %e\n',eps_p,eps_q); 
       end
        
    elseif count == 14

        FLAG1 = strtrim(arg{1});
        
        if strcmp(FLAG1,'SMOOTH_MOD')==0 && strcmp(FLAG1,'SMOOTH_MOD_DIF')==0
            
            fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
            fprintf('Please revise input\n')
            return
        
        end
        fprintf('Gradient: \t\t\t %s\n',FLAG1);
        
    elseif count == 15

        FLAG2 = strtrim(arg{1});
        
        if strcmp(FLAG2,'GRADm')==0 && strcmp(FLAG2,'dmdx')==0
            
            fprintf('FLAG must be GRADm | dmdx\n')
            fprintf('Please revise input\n')
            return
        
        end
        fprintf('Type of derivative: \t\t\t %s\n',FLAG2);
        
    elseif count == 16
        
        temp = strtrim(arg{1});
        
        if isempty(regexp(arg{1},'VALUE','match'))==0
            
            temp = regexp(temp,'\s','split');
            ROT(1) = str2num(temp{2});
            ROT(2) = str2num(temp{3});
            ROT(3) = str2num(temp{4});
        else
            
            ROT = zeros(3,1);
            
        end
        
        fprintf('Objective function rotated azm: %f, dip: %f, plunge: %f\n',ROT(1),ROT(2),ROT(3))
            
    end
    
    line = fgets(fid);
    count = count +1;
    
end
       
fclose(fid);