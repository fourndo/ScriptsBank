function [meshfile,obsfile,m_mag,topofile,wr_flag] = MAG3Csen_read_inp(inputfile)
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
        
        m_mag = strtrim(arg{1});
        if isempty(regexp(m_mag,'null','match'))==0
            m_mag = [];
            fprintf('No magnetization file - assume all induced.\n');
        else
            fprintf('Magnetization model: \t\t %s\n',m_mag);
        end
        
    elseif count == 4
        
        topofile = strtrim(arg{1});
        
        if isempty(regexp(topofile,'null','match'))==0
            topofile = [];
            fprintf('No topography - entire mesh used.\n');
        else
            fprintf('Topography: \t\t %s\n',topofile);
        end
        
    elseif count == 5
        
        wr_flag = strtrim(arg{1});
        fprintf('Weighting: \t\t\t %s\n',wr_flag);
        
    
    end
    
    line = fgets(fid);
    count = count +1;
    
end
       
fclose(fid);