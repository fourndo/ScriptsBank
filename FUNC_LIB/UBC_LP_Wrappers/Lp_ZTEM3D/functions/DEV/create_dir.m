function [file_list,new_dir] = create_dir(p,q,l,home_dir,root_lib)
% Create directory for current inversion parameter. Check if results
% already exist and prompt the option to keep both or overwrite.
            
cd (home_dir);
addpath (root_lib);  

file_list=ls;
new_dir=(['Inv_lp' num2str(p) '_lq' num2str(q) '_lambda' num2str(l)]);
match= [];
for ii=1:size(file_list,1)

    if strcmp(strtrim(file_list(ii,:)),char(new_dir) )==1
        fprintf(['Program has detected existing results for:' new_dir '\n'])
        fprintf('Would you like to overwrite or keep both?\n')
        user_in = input('Press 0 to overwrite, 1 to rename new\n');
        fprintf('\n');
        match=1;
        if user_in==0
            dos (['rmdir  ' new_dir ' /S /Q']);
            dos (['mkdir ' home_dir '\' new_dir]);

        else 
            
            count = 1;
            for jj=1:length(file_list)
             new_dir=(['Inv_lp' num2str(p) '_lq' num2str(q) '_lambda' num2str(l) '_' num2str(count)]);


                if strcmp(strtrim(file_list(ii,:)),char(new_dir) )==1
                    count = count+1;
                    continue
                end

            end
            dos (['mkdir ' home_dir '\' new_dir]);
        end 
   
        
    end

end

if isempty(match)==1
    dos (['mkdir ' home_dir '\' new_dir]);
end

cd (home_dir);