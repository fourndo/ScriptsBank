function data_2_folder(data_dir,work_dir,flag,argin)
% Extract observation files from the data directory and create new
% directory with name specified by the file
% 
% Author: D Fournier
% Last update : April 30th, 2014


%% Driver

home_dir = pwd;

cd(data_dir)

obs_list = ls;

nb_obs = size(obs_list,1)-2;

line_range = argin;

switch line_range 
    case 'all'
        
       range = 1:nb_obs;
       
    otherwise
        
        range = str2num(argin);

end

for ii = range
    

        
        file_name =  strtrim(obs_list(ii+2,:));
    
        
        line_name = regexp(file_name,'(L\d*)(?=_)','match');
        
        if isempty(line_name)==1
            
            continue
            
        end
        line_name = line_name{1};

    
    
    
%     if strcmp(file_name((end-7):(end-4)),'back')==1
        
    dos (['mkdir ' work_dir '\' line_name]);
    
    dos (['copy ' file_name ' ' work_dir '\' line_name '\' file_name '/Y']);
       
    
%     elseif strcmp(file_name((end-10):(end-4)),'forward')==1
%         
%     dos (['mkdir ' output_folder '\' line_name '\Forward']);
%     
%     dos (['copy ' file_name ' ' output_folder '\' line_name '\Forward']);
%     
%     end
     

   
   
end
    
cd(home_dir);
