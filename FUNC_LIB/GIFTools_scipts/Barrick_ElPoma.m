% Barrick processing
clear all
close all

cd C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\ElPoma\GifTools
inp_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\ElPoma\data\Processing';
dsep = '\';

load('proj')
proj = GIFtools(obj);

%% Load mesh
mesh3D(proj); 
mesh = proj.getItem(1); mesh.readFile([inp_dir dsep 'Mesh_15m_Global.msh']);

% Extract nodes
xn = mesh.fetchXnode; xc = mesh.fetchXcenter; nx = length(xc);
yn = mesh.fetchYnode; yc = mesh.fetchYcenter; ny = length(yc);
zn = mesh.fetchZnode; zc = mesh.fetchZcenter; nz = length(zc);

%% Load observation data
MAGinversion.readDobs(proj,[inp_dir dsep 'obs.mag']);
data = proj.getItem(2);
d    = data.data(:,4);
obsx = data.data(:,1);
obsy = data.data(:,2);

%% Load topo
TOPOdata(proj);
topo = proj.getItem(3); topo.readFile([inp_dir dsep 'Barrick_topo.dat']);

%% Create active cell
actv = ACTIVEmodel.createFromTopo(proj,mesh,topo,1); 
actv.writeFile([inp_dir dsep 'nullcell.dat'])

%% Downsample data and extract stats
keep = zeros(size(data)); % Index the cells to keep
sdev = []; % Keep track of standard deviations for group of cells > 5

for jj = 1 : ny
    
    for ii = 1 : nx
        
        % Look for obs inside the cell
        indx = find( obsx >= xn(ii) & obsx <= xn(ii+1) & obsy >= yn(jj) & obsy <= yn(jj+1));
        
        % If found at least one obs inside the cell
        if ~isempty(indx)
            
            % If more than one cell
           if length(indx) > 1
               
               % Find closest cell to center
               [~,id] = min( sqrt( ( obsx(indx) - xc(ii) ).^2 + ( obsy(indx) - yc(jj) ).^2 ) );
               
               keep(indx(id)) = 1;
               
               % If number of obs >=5, then save the standard deviation
               if length(indx) >= 5
                  
                   sdev = [sdev;obsx(indx(id)) obsy(indx(id)) std(d(indx))];
                   
               end
               
           else
               
               keep(indx) = 1;
               
           end
            
        end
        
    end
    
end

figure; scatter(obsx(keep==0),obsy(keep==0),2,'rx'); hold on
scatter(obsx(keep==1),obsy(keep==1),10,'ko'); hold on
scatter(sdev(:,1),sdev(:,2),20,sdev(:,3),'filled');

%% Compute the floor uncertainty from the median of standard devations
floor = median(sdev(:,3));

%% Re-assign elevations from topo
F = scatteredInterpolant(topo.data(:,1),topo.data(:,2),topo.data(:,3));

obs_z = F(obsx(keep==1),obsy(keep==1)) + 2.0;

% Create a new observation item
new_d = MAGdata(proj,{'X','Y','Z','data','Uncert'},...
    [obsx(keep==1) obsy(keep==1) obs_z d(keep==1) ones(sum(keep),1) * floor]);

new_d.setioHeadInd([1:5,0,0]);
new_d.setAinc(data.getAinc);
new_d.setAdec(data.getAdec);
new_d.setInc(data.getAinc);
new_d.setDec(data.getAdec);
new_d.setB0(data.getB0);
new_d.setName('Data_sub_mesh');
MAGinversion.writeDobs(new_d,[inp_dir dsep 'Data_sub_mesh.mag']);
