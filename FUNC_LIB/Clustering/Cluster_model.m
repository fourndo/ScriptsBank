% Load model and compute C-mean clustering
% Dominique Fournier 2015/05/19
% close all
clear all
close all

% addpath ..\FUNC_LIB\;

%% INPUT PARAMETERS
% Project folders
inp_dir = 'Models';

out_dir = 'Cluster_models';

meshfile = 'TKCmesh_sub.msh';


model{1} ='sus_final.sus' ; lbel{1} = 'Susc';
model{2} ='dens_final.den'; lbel{2} = 'Den';

model{3} = 'sig_merged.con';lbel{3} = 'Cond';
% model{3} ='dighemcondtensor.con' ;
% model{4} = 'sig_vtem.con';
model{4} = 'Peta_late.eta' ;lbel{4} = 'Charg1';
model{5} ='Peta_early.eta' ;lbel{5} = 'Charg2';

% Number of cluster
nclusters = 5;

% Convert to LOG10
flag = [0 0 1 0 0];

% Output name
out_name = ['Kmeans_VTEMCon_' num2str(nclusters) 'Clusters.dat'];

% Define cutting planes for viewing model
zsection = 5;

% Define color palette for plotting
cvec = 1:nclusters;
cpal = [255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255;


%% CODE STARTS HERE
% Load mesh file and convert to vectors (UBC format)
[xn,yn,zn] = read_UBC_mesh([inp_dir '\' meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

xc = ( xn(1:end-1) + xn(2:end) ) / 2;
yc = ( yn(1:end-1) + yn(2:end) ) / 2;

mcell = nx*ny*nz;

nmodels = size(model,2);

mmat = zeros(mcell,nmodels);

% Load models
for ii = 1 : nmodels
    
    mmat(:,ii) = load([inp_dir '\' model{ii}]);
    
    % Convert to log space Con model
    if flag(ii) == 1
        
        mmat(:,ii) = log10(mmat(:,ii));
        
    end
    
end

% Remove null cells
nullcell = mmat(:,1) ~= -100;

% Data matrix without nullcell
mmat = mmat(nullcell,:);

mmat_scale = mmat;

% Normalize each column
for ii = 1 : nmodels
    

    
    mmat_scale(:,ii) =( mmat_scale(:,ii) - min(mmat_scale(:,ii))) ;
    mmat_scale(:,ii) =mmat_scale(:,ii) / max(mmat_scale(:,ii));
    
end
% Compute clustering
[idx , c] = kmeans(mmat_scale,nclusters);

% Map back to 3D space
IDX = ones(mcell,1)*-100;

IDX(nullcell) = idx;

% Save result
save([out_dir '\' out_name],'-ascii','IDX');

%% STATS
stat{1} = zeros(nclusters,nmodels);
stat{2} = zeros(nclusters,nmodels);
stat{3} = zeros(nclusters,nmodels);
stat{4} = zeros(nclusters,nmodels);

for ii = 1 : nclusters 
    
    stat{1}(ii,:) = mean(mmat(idx==ii,:),1);
    stat{2}(ii,:) = std(mmat(idx==ii,:),1);
    stat{3}(ii,:) = min(mmat(idx==ii,:),[],1);
    stat{4}(ii,:) = max(mmat(idx==ii,:),[],1);
end

%% Plot model and stats
set(figure, 'Position', [50 25 850 850]); 

m = reshape(IDX,nz,nx,ny);

m2D = reshape(m(zsection,:,:),nx,ny);

ax1 = axes('Position',[0.025 .3 .5 .5]);
h = imagesc(xc,yc,m2D'); hold on
% set(gca,'XTickLabel',[])
set(gca,'YDir','normal')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
axis equal
cc = colormap(cpal(1:nclusters,:));
colorbar('SouthOutside')


for jj = 1 : nmodels
    
    ax1 = axes('Position',[0.5 1 - 0.15*jj .4 .075]);
    
    
for ii = 1 : nclusters
   
    cn = stat{1}(ii,jj);
    ex = stat{2}(ii,jj);
    
    plot([stat{3}(ii,jj) stat{4}(ii,jj)],[ii ii],'k-','LineWidth',2); hold on
    h = fill([cn - ex cn - ex cn + ex cn + ex cn - ex],[ii-0.5 ii+0.5 ii+0.5 ii-0.5 ii-0.5],cpal(ii,:));
    
    text(cn,ii,num2str(ii),'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10)
    
end

    set(ax1,'Box','off','Color','none','YTickLabel',[])
    xlim(ax1 ,[min(mmat(:,jj)) max(mmat(:,jj))]);
    ylim(ax1 ,[0 nclusters+0.5]);
    xlabel(lbel{jj});
    if flag(jj)==1
        
        set(ax1,'xscale','log')
        
    end
end
