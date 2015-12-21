% Function merge tiles of colocated points
clear all
close all

work_dir = 'C:\LC\Private\dominiquef\Projects\4239_Kaminak_Coffee_Mag\Modeling\GRIDDING\Merge';

cd (work_dir)

%% Driver
file_list=ls;
nfile = size(file_list,1);

for jj = 1 : nfile-2
    
    data{jj}=load(strtrim(file_list(jj+2,:)));
    
    % Find center of the tile for weights
    center(1) = mean(data{jj}(:,1));
    center(2) = mean(data{jj}(:,2));
    
    dx = abs(data{jj}(:,1) - center(1)) + 1e-2;
    dy = abs(data{jj}(:,2) - center(2)) + 1e-2;
   
    wx = 1 - dx/max(dx);
    wy = 1 - dy/max(dy);
   % Cosine tapper 
   wx = (-0.5*cos(-wx*pi)+0.5)  ; 
   wy = (-0.5*cos(-wy*pi)+0.5)  ;
   
   w = wx.*wy + 1e-10;
    data{jj} = [data{jj}(:,:) w];
    
    if jj ==1
        cloud = data{jj};
        
    else
        
        cloud = [cloud;data{jj}];
        
    end
end

cloud = [cloud zeros(size(cloud,1),1)];

% Merge all the datasets and average overlapping points
dataout  = zeros(1,4);
count = 1;
progress = -1;
tic
for ii = 1 : size(cloud,1)
    
    % If found new point
    if cloud(ii,6)==0
        
        % Find collocated points
        indexx = cloud(:,1)==cloud(ii,1);
        indexy = cloud(:,2)==cloud(ii,2);
        
        index = (indexx.*indexy)==1;
        
        val = sum(cloud(index,4).*cloud(index,5)) / sum(cloud(index,5));
        
        dataout(count,1:3) = cloud(ii,1:3);
        dataout(count,4)   = val;
        
        % Flag used points
        cloud(index,6) = 1;
        count = count +1;
    end
    
    d_iter = floor(ii/size(cloud,1)*100);
    if  d_iter > progress
        
        fprintf('Computed %i pct of data in %8.5f sec\n',d_iter,toc)
        progress = d_iter;
        
        
         
    end
    
    
end

save([work_dir '\Merged_dataset.xyz'],'-ascii','dataout');