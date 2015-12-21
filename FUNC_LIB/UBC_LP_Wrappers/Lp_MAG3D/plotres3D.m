clear all
close all

% Final results is in np,nq,nl,mcell format
%% Plot in this format:
%           Q
%          ^
%         /
%        /
%       / _ _ _ _ p
%       |
%       |
%       |
%       |
%       lambda (+ up)
%
%% Load results and set up axis vectors
addpath C:\Users\dominiquef\Dropbox\DOM_Projects\INVMAG3D\sliceomatic
work_dir = 'C:\Projects\3796_AGIC_Research\MAG3D\Current\Intrusive\UBC_cluster';
cd (work_dir); cd ..
model = load ('Model_intrusive.sus');
model(model==-1)=-1;

% load MAG3D_model_Intrusive_0p5_v4
% np = length(p);
% nq = length(q);
% nl = length(l);

% load ('C:\Projects\3796_Compactness\MAG3D\Current\Lp_test\Workspace\MAG3D_model_Intrusive.mat');
% Lp-norm on model
q = (0:25:200)/100;

% Lp-norm on gradient (l2-norm by default)
p = (0:25:200)/100;

% Scale phi_cxyy and phi_cm
l = 0:0.25:1.75;

np = length(p);
nq = length(q);
nl = length(l);

cd (work_dir);
folder_list=ls;

range = 1:size(folder_list,1);

finalmodels = zeros(np,nq,nl,length(model));
finalphid = zeros(np,nq,nl,1);

% Cycle through all the files
for oo=range
    
    folder_name = strtrim(folder_list(oo,:));
    if length(folder_name) >=3
    
    
        cd (folder_list(oo,:));
        
    file_list=ls;
    getp =regexp(folder_name,'lp+(.)+lq','tokens');
    ip =  find(p==str2num(char(getp{1})));

    getq =regexp(folder_name,'lq+(.)+lambda','tokens');
    iq =   find(q==str2num(char(getq{1})));
    
    getl =regexp(folder_name,'(?=(a))(.)','split');
    il =   find(l==str2num(char(getl{end})));
    
    finalmodels(ip,iq,il,:) = load (file_list(end,:));
    
% Extract final data misfit from log file
%     fid=fopen('maginv3d.log','rt');

    max_num_lines = 30000;
% Go through the log file and extract data and the last achieved misfit
%         for ii=1:max_num_lines         	
%         line=fgets(fid); %gets next line 
% 
%             if line==-1
%                 fprintf('File ended at line %i\n',ii);
%                 fprintf('Did not find the information needed - review log file\n')
%                 break
%             end
% 
%             if length(strtrim(line))>=length('data misfit:')
%                 description = strtrim(line);
%                 if strcmp(description(1:12),'data misfit:')==1
%                     finalphid(ip,iq,il,:) = str2num(description(13:end));
%                     fclose(fid);
%                     break
%                 end
%             end
%         end
    cd ..
   
    end
    
    
    
end

%% Calculate model norm
L2 = zeros(np,nq,nl);
RMS = zeros(np,nq,nl);
L1 = zeros(np,nq,nl);
Lfro = zeros(np,nq,nl);
Linf = zeros(np,nq,nl);
RMS = zeros(np,nq,nl);
phid = zeros(np,nq,nl);
for ip = 1:np
   for iq = 1:nq
      for il = 1:nl
         rec = finalmodels(ip,iq,il,:);
         rec = rec(:);
         % Norms:
         L2(ip,iq,il) = norm((rec - model),2);
         RMS(ip,iq,il) = sqrt(sum(rec - model).^2/length(model));
         L1(ip,iq,il) = norm((rec - model),1);
         Lfro(ip,iq,il) = norm((rec - model),'fro');
         Linf(ip,iq,il) = norm((rec - model),Inf);
         % Phid:
%          phid(ip,iq,il) = finalphid(iq,ip,il);
      end
   end
end

%% Can also use slicomatic:
sliceomatic(L2,p,q,l);

[p,q,l] = meshgrid(p,q,l);

%% Plot results:
figure
slice(p,q,l,phid,1,1,1)
hold on
slice(p,q,l,phid,1,2,1)
title('\Phi_d');
xlabel('P');
zlabel('\lambda');
ylabel('Q');
axis equal
axis tight
colorbar

figure
slice(p,q,l,L2,1,1,1)
hold on
slice(p,q,l,L2,1,2,1)
caxis([22 25])
title('L-2 norm');
xlabel('P');
zlabel('\lambda');
ylabel('Q');
axis equal
axis tight
colorbar

figure
slice(p,q,l,L1,1,1,1)
caxis([1000 1400])
title('L-1 norm');
xlabel('P');
zlabel('\lambda');
ylabel('Q');
axis equal
axis tight
colorbar

% subplot(2,3,4)
% slice(p,q,l,RMS,1,1,1)
% title('RMS model error');
% xlabel('P');
% zlabel('\lambda');
% ylabel('Q');
% axis equal
% axis tight
% colorbar

% subplot(2,3,5)
% slice(p,q,l,Linf,1,1,1)
% title('L-\infty norm');
% xlabel('P');
% zlabel('\lambda');
% ylabel('Q');
% axis equal
% axis tight
% colorbar

% subplot(2,3,6)
% slice(p,q,l,Lfro,1,1,1)
% title('Frobenius norm');
% xlabel('P');
% zlabel('\lambda');
% ylabel('Q');
% axis equal
% axis tight
% colorbar

