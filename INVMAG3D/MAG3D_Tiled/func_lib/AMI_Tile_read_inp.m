function [meshfile,obsfile,topofile,mstart,mref,m_vec,chi_target,alphas,beta,bounds,lp_MAIvec,lp_MVIvec,lp_tresh,weightfile,FLAG1,FLAG2,tilefile] = AMI_Tile_read_inp(inputfile)
% Function [work_dir,meshfile,obsfile,wr_flag,chi_target,alphas,beta,pvec,qvec,lvec] = MAG3D_read_inp('inputfile')
% Read input file for inversion

bounds = zeros(3,2);
lp_MVIvec = zeros(1,15);
lp_MAIvec = zeros(1,5);
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
            fprintf('Starting model: \t\t %e\n',mref);
            
        else
            
            fprintf('Starting model: \t\t %s\n',mref);
            
        end
         
    elseif count == 6
        
        m_vec = strtrim(arg{1});
        
        if ~isempty(regexp(m_vec,'DEFAULT','match'))
            
            m_vec = [];
            fprintf('Assume inducing direction\n');
            
        else
            
            fprintf('Magnetization vector model: \t\t %s\n',m_vec);
            
        end
        
    elseif count == 7
        
        chi_target = str2num(arg{1});
        
        if isempty(chi_target)==1
            
            chi_target = 1.0;
            
        end
        
        fprintf('Target chi factor: \t %f\n',chi_target);
        
    elseif count == 8
        
        alphas = str2num(arg{1});
        
        if isempty(alphas)==1 || length(alphas)~=4
            
            fprintf('Error in input file at line %i\n',count);
            fprintf('Requires four numerical values (e.g.--> 0.001 1 1 1\n')
            break
            
        end
        
        fprintf('Alpha values: \t\t %4.1e %4.1e %4.1e %4.1e\n',alphas);
        
    elseif count == 9
        
       beta = str2num(arg{1});   
       
       if isempty(beta)==1
            
          fprintf('Starting Beta: \t\t COMPUTED\n'); 
          
       else
           
          fprintf('Starting Beta: \t\t %f\n',beta); 
       end
    
   elseif count == 10
        
        temp = strtrim(arg{1});
        
        if isempty(regexp(temp ,'VALUE','match'))==0
            
            temp = regexp(temp,'\s','split');
            bounds(1,1) = str2num(temp{2});
            bounds(1,2) = str2num(temp{3});
            fprintf('Bounds: [ %f %f] \n',bounds(1,1),bounds(1,2));
            
        else
            
            fprintf('Default bounds on p-component: [-1 -1]\n');
            bounds(1,:) =[-1 1];
            
        end
        
    elseif count == 11
        
        temp = strtrim(arg{1});
        
        if isempty(regexp(temp ,'VALUE','match'))==0
            
            temp = regexp(temp,'\s','split');
            bounds(2,1) = str2num(temp{2});
            bounds(2,2) = str2num(temp{3});
            fprintf('Bounds: [ %f %f] \n',bounds(2,1),bounds(2,2));
            
        else
            
            fprintf('Default bounds on s-component: [-1 -1]\n');
            bounds(2,:) =[-1 1];
            
        end
        
    elseif count == 12
        
        temp = strtrim(arg{1});
        
        if isempty(regexp(temp ,'VALUE','match'))==0
            
            temp = regexp(temp,'\s','split');
            bounds(3,1) = str2num(temp{2});
            bounds(3,2) = str2num(temp{3});
            fprintf('Bounds: [ %f %f] \n',bounds(3,1),bounds(3,2));
            
        else
            
            fprintf('Default bounds on t-component: [-1 -1]\n');
            bounds(3,:) =[-1 1];
            
        end
    
    elseif count == 13
        

        if isempty(regexp(arg{1},'VALUE','match'))==0
            
       	    temp = regexp(arg{1},'\s','split');
            for ii = 2 : 6


                    lp_MAIvec(1,ii-1) = str2num(temp{ii});

                    if isempty(lp_MAIvec(1,ii-1))==1

                        fprintf('Error in input file at line %i\n',count);
                        fprintf('Requires five numerical value (e.g.--> p ,qx, qy, qz , r )\n')
                        break

                    end

            end
            
            fprintf('Lp-norm for the amplitude inversion: %2.1f %2.1f %2.1f %2.1f %2.1f\n',lp_MAIvec);
       
        elseif isempty(regexp(arg{1},'FILE','match'))==0
            
            temp = regexp(arg{1},'\s','split');
            lp_MAIvec{1} = strtrim(temp{2});
            	
        end
        
    elseif count == 14
        

        if isempty(regexp(arg{1},'VALUE','match'))==0
            
       	    temp = regexp(arg{1},'\s','split');
            for ii = 2 : 6


                    lp_MVIvec(1,ii-1) = str2num(temp{ii});

                    if isempty(lp_MVIvec(1,ii-1))==1

                        fprintf('Error in input file at line %i\n',count);
                        fprintf('Requires five numerical value (e.g.--> p ,qx, qy, qz , r )\n')
                        break

                    end

            end
       
            fprintf('Lp-norm for the MVI p-component: %2.1f %2.1f %2.1f %2.1f %2.1f\n',lp_MVIvec(1,1:5));
            
        elseif isempty(regexp(arg{1},'FILE','match'))==0
            
            temp = regexp(arg{1},'\s','split');
            lp_MVIvec{1} = strtrim(temp{2});
            	
        end
        
    elseif count == 15
        
        if isempty(regexp(arg{1},'VALUE','match'))==0
            
       	    temp = regexp(arg{1},'\s','split');
            for ii = 2 : 6


                    lp_MVIvec(1,ii+4) = str2num(temp{ii});

                    if isempty(lp_MVIvec(1,ii+4))==1

                        fprintf('Error in input file at line %i\n',count);
                        fprintf('Requires five numerical value (e.g.--> p ,qx, qy, qz , r )\n')
                        break

                    end

            end
       
            fprintf('Lp-norm for the MVI s-component: %2.1f %2.1f %2.1f %2.1f %2.1f\n',lp_MVIvec(1,6:10));
            
        elseif isempty(regexp(arg{1},'FILE','match'))==0
            
            temp = regexp(arg{1},'\s','split');
            lp_MVIvec{2} = strtrim(temp{2});
            	
        end
        
    elseif count == 16
        
         if isempty(regexp(arg{1},'VALUE','match'))==0
            
       	    temp = regexp(arg{1},'\s','split');
            for ii = 2 : 6


                    lp_MVIvec(1,ii+9) = str2num(temp{ii});

                    if isempty(lp_MVIvec(1,ii+9))==1

                        fprintf('Error in input file at line %i\n',count);
                        fprintf('Requires five numerical value (e.g.--> p ,qx, qy, qz , r )\n')
                        break

                    end

            end
       
            fprintf('Lp-norm for the MVI t-component: %2.1f %2.1f %2.1f %2.1f %2.1f\n',lp_MVIvec(1,11:15));
        elseif isempty(regexp(arg{1},'FILE','match'))==0
            
            temp = regexp(arg{1},'\s','split');
            lp_MVIvec{3} = strtrim(temp{2});
            	
         end
    
     elseif count == 17

        temp = regexp(arg{1},'\s','split');


        if strcmp(temp{1},'EZ')

            lp_tresh{1} = 1;
            lp_tresh{2} = [str2num(temp{2}) str2num(temp{3})];
            fprintf('Effective zero model: \t\t\t %f\n',lp_tresh{2}(1));
            fprintf('Effective zero gradient: \t\t\t %f\n',lp_tresh{2}(2));

        else

            lp_tresh{1} = 0;
            lp_tresh{2} = [];

            fprintf('Effective zero determined by L-curve');

        end

    
    elseif count == 18

        weightfile = strtrim(arg{1});
        
        if strcmp(weightfile,'DEFAULT')
            
            weightfile = [];
            fprintf('No Cell base weight file: \t\t\t %s\n',weightfile);
            
        else
            fprintf('Cell base weight file: \t\t\t %s\n',weightfile);
            
        end
        
        
    elseif count == 19

        FLAG1 = strtrim(arg{1});
        
        if strcmp(FLAG1,'SMOOTH_MOD')==0 && strcmp(FLAG1,'SMOOTH_MOD_DIF')==0
            
            fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
            fprintf('Please revise input\n')
            return
        
        end
        fprintf('Gradient: \t\t\t %s\n',FLAG1);
        
    
    elseif count == 20

        FLAG2 = strtrim(arg{1});
        
        if strcmp(FLAG2,'GRADm')==0 && strcmp(FLAG2,'dmdx')==0
            
            fprintf('FLAG must be GRADm | dmdx\n')
            fprintf('Please revise input\n')
            return
        
        end
        fprintf('Type of derivative: \t\t\t %s\n',FLAG2);
        
    
    elseif count == 21
        
        tilefile = strtrim(arg{1});
        
        if isempty(regexp(tilefile,'NONE','match'))==0
            
            tilefile = [];
            fprintf('Single tile inversion\n');
            
        else
            
            fprintf('Single tile: \t\t %s\n',tilefile);
            
        end
%         tilefile = strtrim(arg{1});
        
    end
    
    line = fgets(fid);
    count = count +1;
    
end
       
fclose(fid);