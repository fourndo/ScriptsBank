close all
clear all


home_dir = 'C:\Users\dominiquef\Dropbox\DOM_Projects\MUON3D';
root_lib = ([home_dir '\functions']);

addpath (root_lib) 
addpath ([home_dir '\data']);

%% USER INPUTS



load kernel_10m 
load Wx_10m 
load Wy_10m 
load Wz_10m 
load Ws_10m
load mnull_10m
% load Zdiscretized

[meshfile]=get_UBCmesh('Mesh_10m_nopad.msh');

dx = meshfile(3,meshfile(3,:)~=0)';
dy = meshfile(4,meshfile(4,:)~=0)';
dz = meshfile(5,meshfile(5,:)~=0)';

nx = meshfile(1,1); %size(X,1);    %number of cell in X
ny = meshfile(1,2); %size(X,2);    %number of cell in Y
nz = meshfile(1,3); %size(X,3);    %number of cell in Z

x0 = meshfile(2,1);
y0 = meshfile(2,2);
z0 = meshfile(2,3);


mcell = nx * ny * nz;

% Make Ws matrix
count = 1;
ws = zeros(mcell,1);

for jj = 1 : ny
    for ii = 1 : nx
        for kk = 1 : nz
            
            ws(count) = dz(kk) * dx(ii) * dy(jj);
            count = count+1;
            
        end
    end
end

% Ws = spdiags(ws,0,mcell,mcell);

model = importdata('Model_Blockcave_nopad.den');
% topo = importdata('Topo_Gaussian.topo');
% surf = topo.data;
% obs_loc = importdata('Obs_loc.obs');


% 2D problem

ndata=size(G,1);

%% Create Noisey data:

% Generate raw data
data= G * model;

pct_noise = 0.05;
noise = (pct_noise.*max(abs(data))).*randn(ndata,1);
% load noise
Wd = spdiags(1./abs(noise),0,ndata,ndata);
d = data + noise;

save('noise','noise')
%% Sensitivity Matrices
% Normalize d and G by standard deviation:
target = sum((Wd*(data - d)).^2);
G = Wd*G;
GtG = G'*G;
d = Wd*d;

phid = @(x) sum((G*x - d).^2);

% Length scale
Lx = 1;%4*min(dx); % Length Scale
Ly = 1;%4*min(dy); % Length Scale
Lz = 1;%4*min(dz); % Length Scale

alphas = 1 / (4*min(dx)).^2 ;
alphax = 1.0;
alphay = 1.0;
alphaz = 1.0;


lvec= 1.8;
pvec= 2;
qvec= 0; 

nl=length(lvec);
nq=length(qvec);
np=length(pvec);
phi_d=99999;
ii = 1;
tic
for pp=1:length(pvec)
    
for qq=1:length(qvec)   
    
for ll=1:length(lvec)
 
    [file_list,new_dir] = create_dir(pvec(pp),qvec(qq),lvec(ll),home_dir,root_lib);
    
    beta = 1e+5;
    invmod = ones(size(model))*1e-4;
    invmod(mnull==0)=0;
    fprintf('Iteration %i of %i.\n',sub2ind([nl,nq,np],ll,qq,pp),np*nq*nl); 
    
    while phi_d(end)>target
   

            if ii==1                 %First iteration

                lambda = 1;

                modfunc= phim_3D(invmod,2,Ws,Wx,Wy,Wz,alphas,alphax,alphay,alphaz,lambda,2,ii);
                
            else
                
                lambda = lvec(ll);

                modfunc= phim_3D(invmod,pvec(pp),Ws,Wx,Wy,Wz,alphas,alphax,alphay,alphaz,lambda,qvec(qq),ii);
            end


    % CG solver
    A=(GtG + beta(end)*(modfunc));
    RHS=(G')*d;
    % [invmod,iter]=CG_solver(invmod,A,RHS);
    invmod=A\RHS;

    phi_d(ii) = phid(invmod)


        if ii>1

            if(abs(phi_d(ii)-phi_d(ii-1))/phi_d(ii) <= 0.01 ) && (ii > 3)
                break;    
            end

        end


        if phi_d(end) < target*2
          beta = 0.75*beta;
        else
          beta = 0.75*beta;
        end
        
    
%         save([new_dir '\model_iter'],'-ascii','invmod')
        ii = ii+1;
        
    end
    
mnull=mnull(:);
invmod(mnull==0)=-100;
phi_d = 99999;
ii = 1;
save(['results/lambda' num2str(lvec(ll)) 'q' num2str(qvec(qq))  'p' num2str(pvec(pp)) '.den'],'-ascii','invmod')
%          finalphid(pp,qq,ll) = phi_d(end);
%          finalmodels(pp,qq,ll,:) = invmod;

end

end

end

fprintf('Results took %6.2f minutes.\n',toc/60);
% save('2Dresults_model_E_scaleALL.mat','finalmodels','lvec','pvec','qvec','finalphid');

