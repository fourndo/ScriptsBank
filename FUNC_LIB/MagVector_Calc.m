% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Modeling\Inversion\Tile_50m_VOI\All_Tiles';
% vec_file = 'Tile31_MVI.fld';
% msh_file = 'Tile31.msh';

tile_file = '..\..\DepositTilesVOI';
pipe_file = '..\..\Pipes43101';
loc_file = '..\..\XYloc.dat';
dsep = '\';

obsfile = '..\..\Obs_Paul_Lake_SUB_1pc_10nT.dat';

dir_ls = ls(work_dir);

%% Load data and plot

[H, avg_I, Dazm, avg_D, Obsx, Obsy, Obsz, data, wd_full] = read_MAG3D_obs([work_dir dsep obsfile]);

xx = min(Obsx):50:max(Obsx);
yy = min(Obsy):50:max(Obsy);
[YY,XX] = ndgrid(yy,xx);

F = scatteredInterpolant(Obsy, Obsx, data ,'natural');

grid_d = F(YY,XX);

h = set(figure, 'Position', [50 25 1500 1000]);
% msh = mesh(Xn,Yn,ones(size(Xn)),'FaceColor','none'); hold on
view([0 90]) 
imagesc(xx,yy,grid_d);hold on
caxis([-500 250])
colormap(jet)
xlim([min(Obsx) max(Obsx)]); rangex = max(Obsx) - min(Obsx); 
ylim([min(Obsy) max(Obsy)]); rangey = max(Obsy) - min(Obsy);
set(gca,'Ydir','normal')
axis equal tight

%% Load entire model
[xn,yn,zn] = read_UBC_mesh([work_dir dsep '..\..\Mesh_50m_v2.msh']);
xc = ( xn(2:end) + xn(1:end-1) )/2; nx = length(xc);
yc = ( yn(2:end) + yn(1:end-1) )/2; ny = length(yc);
zc = ( zn(1:end-1) + zn(2:end) )/2; nz = length(zc);
mcell = nx*ny*nz;

[Xc,Yc] = ndgrid(xc,yc);

mvec = load([work_dir dsep 'Merged_M_model.fld']);
mamp = load([work_dir dsep 'Merged_Amp_model.amp']);
mamp = reshape(mamp,nz,nx,ny);

nullcell = mamp~=-100;

%% Load tiles and deposit name
load([work_dir dsep tile_file]); tiles = DepositTilesVOI;
load([work_dir dsep pipe_file]); pipes = Pipes43101;

load([work_dir dsep loc_file])

fid = fopen([work_dir dsep '..\Pipes_M_Stats.dat'],'w');
fprintf(fid,'%13s &\tEasting (m) &\tNorthing (m) &\t $\\kappa_{e}$  &\t Incl (dd.d) &\t Decl (dd.d) \\\\ \n','Pipe_ID');
 
count = 0;
for ii = 1 : size(tiles,1)


if XYloc(ii,3) ~= 2
    
    count = count + 1;
    
else
    
    continue
    
end

% [xn,yn,zn] = read_UBC_mesh([work_dir dsep 'Tile' num2str(count) '.msh']);
% xc = ( xn(2:end) + xn(1:end-1) )/2; nx = length(xc);
% yc = ( yn(2:end) + yn(1:end-1) )/2; ny = length(yc);
% zc = ( zn(1:end-1) + zn(2:end) )/2; nz = length(zc);
% mcell = nx*ny*nz;
% 
% [Xc,Yc] = ndgrid(xc,yc);
% 
% mvec = load([work_dir dsep 'Tile' num2str(count) '_MVI.fld']);
% mamp = load([work_dir dsep 'Tile' num2str(count) '_MAI_esus.sus']);
% mamp = reshape(mamp,nz,nx,ny);
% 
% nullcell = mamp~=-100;



% M = sqrt(sum(mvec.^2,2));

% Find deposit in the database
for jj = 1 : size(pipes,1)
    
    if strcmp(pipes{jj,1},tiles{ii,5})
        
        fprintf('Pipe found: %s\n',pipes{jj,1})
        
        break
        
    end
    
end

% % Select the cells closest XY to deposit location
% R = sqrt( (Xc - pipes{jj,5}).^2 + (Yc - pipes{jj,6}).^2 );
% 
% [~,id2d] = min(R(:));
% [idx,idy] = ind2sub([nx ny],id2d);
% 
% % Find the first cell below topo
% idz = find( mamp(:,idx,idy) ~= -100 , 1);

%% Load section through the amplitude model and select true center
% if ii ==55
% [~, ~, ~, ~, obsx, obsy, obsz, d, ~] = read_MAG3D_obs([work_dir dsep 'Tile' num2str(ii) '_MVI.pre']);
% slice2D = permute( squeeze(mamp(idz,:,:)), [2 1]);
% 
% h = set(figure(2), 'Position', [50 25 1500 1000]);
% subplot(2,1,1)
% mesh(xc,yc,slice2D,'FaceColor','interp',...
%    'EdgeColor','none');hold on
% view([0 90])
% caxis([0 max(slice2D(:))])
% axis([min(obsx) max(obsx) min(obsy) max(obsy)])
% axis square
% title(pipes{jj,1})
% x = min(obsx):50:max(obsx);
% y = min(obsy):50:max(obsy);
% 
% [X,Y] = meshgrid(x,y);
% 
% d_interp     = griddata(obsx, obsy, d,X,Y,'linear'); 
% % data_interp(isnan(data_interp)) = min(data_interp(:))*2;
% 
% 
% subplot(2,1,2)
% h =imagesc(x,y,d_interp);hold on
% set(h,'alphadata',~isnan(d_interp))
% caxis([min(d_interp(:)) max(d_interp(:))]);
% colormap(jet);
% scatter(obsx,obsy,2,'k.')
% set(gca,'YDir','normal')
% axis([min(obsx) max(obsx) min(obsy) max(obsy)])
% axis square
% 
% % Re-define the pipe location with ginput
% idxy(ii,:) = ginput(1);
% keeper(ii) = waitforbuttonpress;
% fprintf('wow\n');
% end
% load([work_dir '..\..\..\XYloc.dat']);
%% Grab a block of cells around the approximate location

[~,idx] = min( abs(xc - XYloc(ii,1)));
[~,idy] = min( abs(yc - XYloc(ii,2)));
% Find the first cell below topo
idz = find( mamp(:,idx,idy) ~= -100 , 1);

if idx+1 > nx || idy+1>ny
    
    continue
    
end
[idz,idx,idy] = ndgrid(idz:idz+2,idx-2:idx+2,idy-2:idy+2);

indx = sub2ind([nz nx ny],idz(:),idx(:),idy(:));

mvec_sub = mvec(indx,:);

% Remove zeros
indx = sum(mvec_sub,2) ~= 0;
mvec_sub = mvec_sub(indx,:);

M =  sqrt(sum(mvec_sub.^2,2));

if size(mvec_sub,1) > 1
% Compute stats on recovered vector within VOI
avg_m = sum(mvec_sub.*repmat(M,[1,3]))/sum(M);
std_m = std(mvec_sub);

end

M = sqrt(sum(avg_m.^2));
std_M = sqrt(sum(std_m.^2));

if M == 0
  
    continue
    
end
% normvec = 1./M;

mx = avg_m(1);
my = avg_m(2);
mz = avg_m(3);

% Normalize the mag vector
mx = mx / M;
my = my / M;
mz = mz / M;

% Unit vector
% for ii = 1 : mcell

avg_I = asind(mz);
avg_D = atand(my./mx);

std_I = std( asind( mvec_sub(:,3)./sqrt(sum(mvec_sub.^2,2) ) ) );
Decl = atand(mvec_sub(:,2)./mvec_sub(:,1) ) ;

% Modify azimuth from North to Cartesian
    
avg_D(my < 0 & mx<0) = -avg_D(my < 0 & mx<0)-90;
avg_D(my > 0 & mx<0) = avg_D(my > 0 & mx<0)-180;
avg_D = mod(450-avg_D,360);

Decl(mvec_sub(:,2) < 0 & mvec_sub(:,1)<0) = -Decl(mvec_sub(:,2) < 0 & mvec_sub(:,1)<0)-90;
Decl(mvec_sub(:,2) > 0 & mvec_sub(:,1)<0) = Decl(mvec_sub(:,2) > 0 & mvec_sub(:,1)<0)-180;
Decl = mod(450-Decl,360);

% Extra change for 360|0 region
if ( sum(Decl > 270 & Decl < 360) * sum(Decl < 90) ) ~=0
    
    Decl(Decl<90) = Decl(Decl<90)+360;
    
end
std_D = std(Decl);

% Normalize standard deviation by the projected length of the XY vector
std_D = std_D * (sqrt(sum(avg_m(1:2).^2)) / M);

%% Plot location and info about M
m_pow = round(log10(M));
m_exp = M / 10^m_pow;

figure(1)


% Convert magnitude to exponential form

if XYloc(ii,3) ~= 2
    
    plot(XYloc(ii,1),XYloc(ii,2),'ro');
    text(XYloc(ii,1),XYloc(ii,2),[pipes{jj,1}],'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(XYloc(ii,1),XYloc(ii,2),['M: ' num2str(round(m_exp*10)/10) 'e' num2str(m_pow)],'VerticalAlignment','top','HorizontalAlignment','left')
    text(XYloc(ii,1),XYloc(ii,2)-400,['I: ' num2str(round(avg_I*10)/10)],'VerticalAlignment','top','HorizontalAlignment','left')
    text(XYloc(ii,1),XYloc(ii,2)-800,['D: ' num2str(round(avg_D*10)/10)],'VerticalAlignment','top','HorizontalAlignment','left')
    
    % Write to file
    fprintf(fid,'%13s &\t%6.0f &\t%7.0f &\t%6.1e $\\pm$ %6.1e &\t%6.1f $\\pm$ %6.1f &\t%6.1f $\\pm$ %6.1f \\\\ \n',...
    pipes{jj,1},XYloc(ii,1),XYloc(ii,2),...
    M,std_M,avg_I,std_I,avg_D,std_D);

else
    
    plot(XYloc(ii,1),XYloc(ii,2),'rx');
%     text(XYloc(ii,1),XYloc(ii,2),[pipes{jj,1}],'VerticalAlignment','bottom','HorizontalAlignment','left')
end


    
% annotation('textbox',...
%     [(pipes{jj,5}-min(Obsx))/rangex (pipes{jj,6}-min(Obsy))/rangey 0.1 0.1],...
%     'String',{pipes{jj,1},...
%     ['M: ' num2str(round(m_exp*10)/10) 'e' num2str(m_pow)],...
%     ['I: ' num2str(round(I*10)/10)],...
%     ['D: ' num2str(round(D*10)/10)]},...
%     'FontSize',10)
    
end

fclose(fid)