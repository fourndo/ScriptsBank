% Tile_Builder
% 
% Function takes in a mesh and data points and designed a tiling scheme
% with constraints on the overlap and data coverage.
% Used for the MAG3D_Tile code.
%
% INPUT
% meshfile
% Obsfile
% Olap
% 
% OUTPUT
% Tile file: SW-NE corner of each tile [x1 y1 x2 y2; ...]

clear all
close all

%% INPUT VARIABLES
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion';

meshfile = 'Mesh_50m.msh';

obsfile = 'Obs_Paul_Lake_SUB_2pc_10nT_DETREND.dat';

dsep = '\';

max_mcell = 1e+5;
min_Olap = 5e+2;
%% Load in files
[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);
dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);

dx = min(dx);
dy = min(dy);

[Xn,Yn] = ndgrid(xn,yn);


%% Plot mesh and obs

h = set(figure, 'Position', [50 25 1500 1000]);
[H, I, Dazm, D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);

load([work_dir dsep 'Obs_flag']);

F = scatteredInterpolant(Obsy, Obsx, data ,'natural');

grid_d = F(YY,XX);

h = set(figure, 'Position', [50 25 1500 1000]);
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
view([0 90]) 
h = imagesc(xx,yy,grid_d);hold on
set(h,'alphadata',~isnan(grid_d) .* flag)
caxis([-250 250])
colormap(jet)
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight

hold on

xlabel('Easting (m)')
ylabel('Northing (m)')



%% Load extra data
load([work_dir '\Pipes43101']);
% 
% for ii = 1 : size(Pipes43101,1)
%     
%     plot(Pipes43101{ii,5},Pipes43101{ii,6},'rx')
%     text(Pipes43101{ii,5},Pipes43101{ii,6},Pipes43101{ii,1});
%     
% end
% 
% load([work_dir '\TargetVOI'])
% for ii = 1 : size(TargetVOI,1)
%     
%     plot(TargetVOI{ii,2},TargetVOI{ii,3},'ro','MarkerSize',10)
%     text(TargetVOI{ii,2},TargetVOI{ii,3},TargetVOI{ii,1},'FontSize',15);
%     
% end







%% Begin tiling algorithm in 1D
% First put a tile on both corners of the mesh
% Add intermediate tiles until the min_Olap is respected
% Repeat for x and y-axis

lx = floor(sqrt(max_mcell/nz));

% In the x-direction
ntile = 1;
Olap  = -1;

while Olap < min_Olap
    
    ntile = ntile + 1;
    % Set location of SW corners
    x0 = [xx(1) xx(end) - lx*dx];
    
    dx_t = round( ( x0(2) - x0(1) ) / ( (ntile-1) *dx) ) ;

    x1 = [x0(1) xx(1) + cumsum((ones(1,ntile-2)) * dx_t * dx) x0(2)];
    x2 = x1 + lx*dx;

    y1 = ones(1,length(x1))*yy(1);
    y2 = ones(1,length(x1))*(yy(1) + lx*dx);
    
    Olap = x1(1) + lx*dx - x1(2);

end
% patch([x1;x1;x2;x2],[y1;y2;y2;y1],[0.5 0 0],'LineWidth',2)
% alpha(0.25);

% Save x-corner location
xtile = [x1' x2'];

% In the y-direction
ntile = 1;
Olap  = -1;

ly = floor(sqrt(max_mcell/nz));

while Olap < min_Olap
    
    ntile = ntile + 1;
    % Set location of SW corners
    y0 = [yy(1) yy(end) - ly*dy];
    
    dy_t = round( ( y0(2) - y0(1) ) / ( (ntile-1) *dy) ) ;

    y1 = [y0(1) yy(1) + cumsum((ones(1,ntile-2)) * dy_t * dy) y0(2)];
    y2 = y1 + ly*dy;

    x1 = ones(1,length(y1))*xx(1);
    x2 = ones(1,length(y1))*(xx(1) + ly*dy);
    
    Olap = y1(1) + ly*dy - y1(2);

end


% Save x-corner location
ytile = [y1' y2'];

%% Replicate the tiling in 2D
Xtile = kron( ones(size(ytile,1),1) , xtile );
Ytile = kron( ytile , ones(size(xtile,1),1) );

%% Generate tiles from deposit location
% tiles = load([work_dir dsep 'XYloc.dat']);
% 
% lx = floor(sqrt(max_mcell/nz)/2);
% ntile = 0;
% Xtile = [];
% Ytile = [];
% 
% for ii = 1 : size(tiles,1)
%         
%     ntile = ntile + 1;
%     
%     xx = tiles(ii,1);
%     yy = tiles(ii,2);
%     
%     % Set location of SW corners
%     x0 = [xn(1) xn(end) - lx*dx];
%     
%     dx_t = round( ( x0(2) - x0(1) ) / ( (ntile-1) *dx) ) ;
% 
%     Xtile = [Xtile ;xx - lx * dx xx + lx * dx];
% 
%     Ytile = [Ytile ; yy - lx * dy yy + lx * dy];
% 
% end

%% Remove redundent tiles

kill = ones(size(Xtile,1),1);

for ii = 1 : size(Xtile,1)
    
    idx = Obsx > Xtile(ii,1) & Obsx < Xtile(ii,2) & Obsy > Ytile(ii,1) & Obsy < Ytile(ii,2);
    
    ndata = sum(idx);
    
    fprintf('Tile %i, ndata: %i\n',ii,ndata);
    if ndata == 0 %|| tiles(ii,3) == 0
        
        kill(ii) = 0;
        
    end
    
    
end

Xtile = Xtile(kill==1,:);
Ytile = Ytile(kill==1,:);

% Plot tiles
h=patch([Xtile(:,1)';Xtile(:,1)';Xtile(:,2)';Xtile(:,2)'],[Ytile(:,1)';Ytile(:,2)';Ytile(:,2)';Ytile(:,1)'],[0.5 0 0],'LineWidth',1);
alpha(h,0);



%% Allow to add tiles manually
xc = 1;
while ~isempty(xc)
    
    [xc,yc] = ginput(1);

    Xnew = [(xc - lx/2*dx) (xc + lx/2*dx)];
    Ynew = [(yc - lx/2*dy) (yc + lx/2*dy)];

    if isempty(xc)

        
        continue

    end
    patch([Xnew(1);Xnew(1);Xnew(2);Xnew(2)],[Ynew(1);Ynew(2);Ynew(2);Ynew(1)],[0 0.5 0],'LineWidth',2,'EdgeColor','r')

    Xtile = [Xtile;Xnew];
    Ytile = [Ytile;Ynew];

end

% Save to tile file
Tiles = [Xtile(:,1) Ytile(:,1) Xtile(:,2) Ytile(:,2)];
save([work_dir dsep 'Tiles_RENAME.dat'],'-ascii','Tiles');

% Tiles = load([work_dir dsep 'Tiles_VOI.dat']);
% patch([Tiles(:,1)';Tiles(:,1)';Tiles(:,3)';Tiles(:,3)'],[Tiles(:,2)';Tiles(:,4)';Tiles(:,4)';Tiles(:,2)'],[0.5 0 0],'LineWidth',1)


%% Plot tile number
for ii = 1: size(Tiles,1)
    
    text(mean(Tiles(ii,[1 3])),mean(Tiles(ii,[2 4])),num2str(ii),'VerticalAlignment','top','HorizontalAlignment','center');
    
end

%% Plot deposit id
% load([work_dir '\DepositTilesVOI']);
% 
% temp = find(kill==1);
% fid = fopen([work_dir '\Deposit_Tiles.dat'],'w');
% for ii = 1 : size(Tiles,1)
%     
%     text(mean(Tiles(ii,[1 3])),mean(Tiles(ii,[2 4])),DepositTilesVOI{temp(ii),5},'VerticalAlignment','bottom','HorizontalAlignment','center');
%     fprintf(fid,'%8.5e\t%8.5e\t%8.5e\t%8.5e\t%s\n',Tiles(ii,1),Tiles(ii,2),Tiles(ii,3),Tiles(ii,4),DepositTilesVOI{temp(ii),5});
%     
% %     text(mean(Tiles(ii,[1 3])),mean(Tiles(ii,[2 4])),Pipes43101{temp(ii),1},'VerticalAlignment','bottom','HorizontalAlignment','center');
% %     fprintf(fid,'%8.5e\t%8.5e\t%8.5e\t%8.5e\t%s\n',Tiles(ii,1),Tiles(ii,2),Tiles(ii,3),Tiles(ii,4),Pipes43101{temp(ii),1});
% end
% 
% fclose(fid);

%% Add colorbar
ax = axes('Position',[0.8 .5 .1 .3]);

% bb = interp1([cvec'],[255 255 255;77 190 238;255 255 0;255 127 0;255 0 0;255 153 200]/255,0:1e-4:mmax);

cbar = colorbar(ax,'EastOutside');
colormap(ax,jet);
caxis([-250 250])
% set(cbar,'Ticks',[0 0.2 0.4 0.6 0.8 1])
% set(cbar,'TickLabels',round(cvec*10000)/10000)

set(gca,'Visible','off');
text(1.5,-.1,'$(nT)$', 'interpreter', 'latex','FontSize',10,'HorizontalAlignment','center')
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
% view([0 90]) 
% scatter(Obsx,Obsy,5,data);hold on
% caxis([-250 250])
% 
% xlim([min(Obsx) max(Obsx)])
% ylim([min(Obsy) max(Obsy)]);
% 
% axis equal