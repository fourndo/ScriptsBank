function [data_sort,xyz_sort,std_sort] = sort_EM1DTM(data,xyz,limits,filter_r,argumin)
% [data_out] = sort_EM1DTM(work_dir,rawdata.mat,radius)
%
% INPUT:
% rawdata:  Matlab array for data
% radius:   Minimum distance between points. Used to filter out data.
%
% OUTPUT:
% data:     Cell array of data with format
% 
% Written by: D.Fournier
% Last update: 2014-03-21

%% FOR DEV ONLY
% clear all
% close all
% 
% work_dir= 'C:\Users\dominiquef\Dropbox\DIGHEM\Codes\Test';
% freqin = [56000 7200 900 5000 900];
% rawdata= 'DIGHEM_data';
% radius = 100;

%% SCRIPT STARTS HERE
% Total number of data
ndat = size(data,1);
ntc = size(data,2);

X = xyz(:,1);
Y = xyz(:,2);

mask = zeros(ndat,1);

if isempty(limits)==1
    
    mask(:) = 1;
    
else
    
    xmax = limits(2);
    xmin = limits(1);

    ymax = limits(4);
    ymin = limits(3);

    mask(xyz(:,1) > xmin & xyz(:,1) < xmax &...
        xyz(:,2) > ymin & xyz(:,2) < ymax) = 1 ;

end
figure; scatter(X,Y);title('Before sorting');hold on


data_sort = zeros(ndat,ntc);
std_sort = zeros(ndat,ntc);
xyz_sort = zeros(ndat,4);

count = 1;
for ii = 1 : ndat

    if mask(ii)==1
%         temp = zeros(length(rc_specs),1); temp(ii:ii+nfreq-1)=1;
%             temp = stn_num(:,1)==stn(ii,1);
        r = ( (X(ii) - X(:)).^2 +...
            (Y(ii) - Y(:)).^2 ) .^0.5;

        % Only keep the curretn points and neighbours at distance r+
        mask(r <= filter_r) = 0;
        mask(ii) = 1;
% 
%         % Select data around and stack
%         index = r < interp_r;
% 
%         % Extract data around location
%         stck = data(index,:);
% 
%         % Pre-allocate matrix
%         avgV = zeros(1,ntc);
% 
%         % Computed weighted average
%         for jj = 1 : ntc
% 
%              % Only keep positive data for each time channel
%             %          stck(:,jj) = abs(stck(:,jj));
%              select = stck(:,jj)>0;
% 
%              std_sort(count,jj) = std(stck(:,jj));
% 
%              mu = mean(stck(:,jj));
%              w = 1./abs(stck(select,jj) - mu );
%              avgV(jj) = sum(w.*stck(select,jj)) / sum(w);
% 
% 
% %              semilogx(tc(jj),mu,'ro'); hold on
% %              semilogx(tc(jj),stck(:,jj),'*','MarkerSize',2); hold on
% 
%         end
% 
%         data_sort(count,:) = avgV;
%         xyz_sort(count,:) = xyz(ii,:);
%         
%         count = count+1;
        
    end

end

data_sort   = data(mask==1,:);
xyz_sort    = xyz(mask==1,:);
std_sort    = abs(data_sort)*0.1;

switch argumin
    
    case 'pos'
 
        tester = data_sort>0;
        
        tester = sum(tester,2);
        
%         mask = tester < size(data,2)/2;
        mask = tester < 2;
        
        data_sort = data_sort(mask==0,:);
        xyz_sort = xyz_sort(mask==0,:);
        std_sort = std_sort(mask==0,:);
        
        fprintf('%i data points were removed because of large negative value\n', sum(mask==1));
            
        
end

figure(1)
scatter(xyz_sort(:,1),xyz_sort(:,2),'r*');title('After sorting')