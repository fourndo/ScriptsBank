% Corrupt 3D data with Gaussian noise
clear all
close all

data = importdata('mutofwd3d.pre');

ndata = size(data.data,1); 
pct_noise = 0.1;
noise = (pct_noise.*max(abs(data.data(:,8)))).*randn(ndata,1);

d = data.data(:,8) + noise;

% Build new file with corrupted data and noise
data_noisy=[data.data(:,1:7) d noise];

%% Write output
fid = fopen('data_10pct_noise.dat','wt');
fprintf(fid,'%i\n',ndata);
for ii=1:ndata
   fprintf(fid,'%i %10.3f %10.3f %10.3f %10.3f %10.3f  %10.3f  %10.6e  %10.6e\n',...
       data_noisy(ii,1),data_noisy(ii,2),data_noisy(ii,3),data_noisy(ii,4),...
       data_noisy(ii,5),data_noisy(ii,6),data_noisy(ii,7),data_noisy(ii,8),...
       data_noisy(ii,9));
end
fclose(fid);
