% Create figure through 3D model
clear all
close all


addpath ..\..\FUNC_LIB;

%% Input Files
work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\TKC\DIGHEM_TMI\25m';

meshfile = 'Mesh_25m.msh';

model2 = 'MAG3D_TMI_iter_7.sus';
model1 = 'MAG3D_TMI_iter_2.sus';

norm_vec = 'Lp_vec.txt';
wmodel = 'Wvec.txt';

zpanel = 5;
ypanel = 10;

padE = 6;
padW = 6;

padN = 6;
padS = 6;

padT = 0;
padB = 12;

iso_cut = 0.00175;
obsfile = '..\DIGHEM_Mag_2pc_floor10nt_25m_ROT.obs';
predfile_l2 = 'MAG3D_TMI_iter_2.pre';
predfile_lp = 'MAG3D_TMI_iter_7.pre';

cam_ang = [30 60];
%% Load in model and plot
set(figure, 'Position', [50 0 900 1000]); 

% Load data
[H, I, Dazm, D, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);


[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

% Move to local coordinates
obsx = obsx - xn(1);
obsy = obsy - yn(1);
% obsz = obsz - zn(1);

xn = xn - xn(1);
yn = yn - yn(1);
% zn = zn - zn(1);

dx = xn(2:end) - xn(1:end-1); nx = length(dx);
dy = yn(2:end) - yn(1:end-1); ny = length(dy);
dz = zn(1:end-1) - zn(2:end); nz = length(dz);
xx = (xn(2:end) + xn(1:end-1))/2;   xx = xx(padW+1:end-padE);
yy = (yn(2:end) + yn(1:end-1))/2;   yy = yy(padS+1:end-padN);
zz = (zn(2:end) + zn(1:end-1))/2;   zz = zz(padT+1:end-padB);

[XX,ZZ,YY] = meshgrid(xx,zz,yy); 
%%
m1 = load([work_dir '\' model1]);
m1 = reshape(m1,nz,nx,ny);

% Get top of model
m2D = get_model_top(m1,nx,ny,nz,-100,1);

temp = reshape(m2D(padW+1:end-padE,padS+1:end-padN),nx-padE-padW,ny-padN-padS); 
temp = (temp');

ax = axes('Position',[0.095 .612 .38 .38]);
h = imagesc(xx,yy,temp);
caxis([0 0.015]);
cmap =[1 1 1;0.900280117988586 0.963585436344147 0.990476191043854;0.800560235977173 0.927170872688293 0.980952382087708;0.700840353965759 0.89075630903244 0.971428573131561;0.601120471954346 0.854341745376587 0.961904764175415;0.501400589942932 0.817927181720734 0.952380955219269;0.401680678129196 0.78151261806488 0.942857146263123;0.301960796117783 0.745098054409027 0.933333337306976;0.328808456659317 0.754902005195618 0.897435903549194;0.355656117200851 0.764705896377563 0.861538469791412;0.382503777742386 0.774509787559509 0.82564103603363;0.40935143828392 0.7843137383461 0.789743602275848;0.436199098825455 0.79411768913269 0.753846168518066;0.463046759366989 0.803921580314636 0.717948734760284;0.489894419908524 0.813725471496582 0.682051301002502;0.51674211025238 0.823529422283173 0.64615386724472;0.543589770793915 0.833333373069763 0.610256433486938;0.570437431335449 0.843137264251709 0.574358999729156;0.597285091876984 0.852941155433655 0.538461565971375;0.624132752418518 0.862745106220245 0.502564132213593;0.650980412960052 0.872549057006836 0.466666668653488;0.677828073501587 0.882352948188782 0.430769234895706;0.704675734043121 0.892156839370728 0.394871801137924;0.731523394584656 0.901960790157318 0.358974367380142;0.75837105512619 0.911764740943909 0.32307693362236;0.785218715667725 0.921568632125854 0.287179499864578;0.812066376209259 0.9313725233078 0.251282066106796;0.838914036750793 0.941176474094391 0.215384617447853;0.865761697292328 0.950980424880981 0.179487183690071;0.892609357833862 0.960784316062927 0.143589749932289;0.919457018375397 0.970588207244873 0.107692308723927;0.946304678916931 0.980392158031464 0.0717948749661446;0.973152339458466 0.990196108818054 0.0358974374830723;1 1 0;1 0.966666638851166 0;1 0.933333337306976 0;1 0.899999976158142 0;1 0.866666674613953 0;1 0.833333313465118 0;1 0.800000011920929 0;1 0.766666650772095 0;1 0.733333349227905 0;1 0.699999988079071 0;1 0.666666686534882 0;1 0.633333325386047 0;1 0.600000023841858 0;1 0.566666662693024 0;1 0.533333361148834 0;1 0.5 0;1 0.466666668653488 0;1 0.433333337306976 0;1 0.400000005960464 0;1 0.366666674613953 0;1 0.333333343267441 0;1 0.300000011920929 0;1 0.266666680574417 0;1 0.233333334326744 0;1 0.200000002980232 0;1 0.16666667163372 0;1 0.133333340287209 0;1 0.100000001490116 0;1 0.0666666701436043 0;1 0.0333333350718021 0;1 0 0];
colormap(cmap);
axis square
hold on
grid on
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
set(gca,'YDir','normal')
hold on
set(gca,'XTickLabel',[])
set(get(gca,'YLabel'),'Rotation',360);
scatter(obsx,obsy,1,'k')
set(gca, 'YAxisLocation', 'right')
text(min(xx)-dx(1), mean(yn),'$(a)$','interpreter', 'latex','FontSize',14);

% Add DO27 marker
scatter3(2000,3035,450,50,'*','r')
text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12);

scatter3(2350,3600,450,50,'*','r')
text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12);

%% Add color bar
%%
ax = axes('Position',[0.095 .215 .38 .38]);
colormap(ax,cmap);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',round([min(temp(:)) median(temp(:)) max(temp(:))]*100)/100)
set(gca,'Visible','off');
text(0.31,1.25,'$Susceptibility\;(SI)$', 'interpreter', 'latex','FontSize',14)

%%
axes('Position',[0.58 .625 .38 .38])
temp =m1( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
aa = isosurface(XX,YY,ZZ,temp,iso_cut);

% Create color map for depth
bb = interp1([max(zz);median(zz);min(zz)],[215 48 39;255 255 191;69 117 180]/255,aa.vertices(:,3));

view(cam_ang);
p = patch('Faces',aa.faces,'Vertices',aa.vertices,'FaceColor','b','LineStyle','none','FaceColor', 'flat','FaceVertexCData',bb);
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
hold on
axis equal
set(gca,'DataAspectRatio',[1 1 1/2])
camlight HEADLIGHT 
lighting phong
set(gca,'YDir','normal')
grid on
xlim([min(xx) max(xx)])
ylim([min(yy) max(yy)])
scatter3(obsx,obsy,obsz,1,'k')

% Add DO27 marker
scatter3(2000,3035,450,50,'*','r')

scatter3(2350,3600,450,50,'*','r')


%%
m2 = load([work_dir '\' model2]);
m2 = reshape(m2,nz,nx,ny);

% Get top of model
m2D = get_model_top(m2,nx,ny,nz,-100);

temp = reshape(m2D(padW+1:end-padE,padS+1:end-padN),nx-padE-padW,ny-padN-padS); 
temp = (temp');


axes('Position',[0.095 .165 .38 .38])
h = imagesc(xx,yy,temp);
caxis([0 0.015]);
cmap =[1 1 1;0.900280117988586 0.963585436344147 0.990476191043854;0.800560235977173 0.927170872688293 0.980952382087708;0.700840353965759 0.89075630903244 0.971428573131561;0.601120471954346 0.854341745376587 0.961904764175415;0.501400589942932 0.817927181720734 0.952380955219269;0.401680678129196 0.78151261806488 0.942857146263123;0.301960796117783 0.745098054409027 0.933333337306976;0.328808456659317 0.754902005195618 0.897435903549194;0.355656117200851 0.764705896377563 0.861538469791412;0.382503777742386 0.774509787559509 0.82564103603363;0.40935143828392 0.7843137383461 0.789743602275848;0.436199098825455 0.79411768913269 0.753846168518066;0.463046759366989 0.803921580314636 0.717948734760284;0.489894419908524 0.813725471496582 0.682051301002502;0.51674211025238 0.823529422283173 0.64615386724472;0.543589770793915 0.833333373069763 0.610256433486938;0.570437431335449 0.843137264251709 0.574358999729156;0.597285091876984 0.852941155433655 0.538461565971375;0.624132752418518 0.862745106220245 0.502564132213593;0.650980412960052 0.872549057006836 0.466666668653488;0.677828073501587 0.882352948188782 0.430769234895706;0.704675734043121 0.892156839370728 0.394871801137924;0.731523394584656 0.901960790157318 0.358974367380142;0.75837105512619 0.911764740943909 0.32307693362236;0.785218715667725 0.921568632125854 0.287179499864578;0.812066376209259 0.9313725233078 0.251282066106796;0.838914036750793 0.941176474094391 0.215384617447853;0.865761697292328 0.950980424880981 0.179487183690071;0.892609357833862 0.960784316062927 0.143589749932289;0.919457018375397 0.970588207244873 0.107692308723927;0.946304678916931 0.980392158031464 0.0717948749661446;0.973152339458466 0.990196108818054 0.0358974374830723;1 1 0;1 0.966666638851166 0;1 0.933333337306976 0;1 0.899999976158142 0;1 0.866666674613953 0;1 0.833333313465118 0;1 0.800000011920929 0;1 0.766666650772095 0;1 0.733333349227905 0;1 0.699999988079071 0;1 0.666666686534882 0;1 0.633333325386047 0;1 0.600000023841858 0;1 0.566666662693024 0;1 0.533333361148834 0;1 0.5 0;1 0.466666668653488 0;1 0.433333337306976 0;1 0.400000005960464 0;1 0.366666674613953 0;1 0.333333343267441 0;1 0.300000011920929 0;1 0.266666680574417 0;1 0.233333334326744 0;1 0.200000002980232 0;1 0.16666667163372 0;1 0.133333340287209 0;1 0.100000001490116 0;1 0.0666666701436043 0;1 0.0333333350718021 0;1 0 0];
colormap(cmap);
% colorbar('NorthOutside')
axis square
hold on
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
hold on
grid on
set(gca,'YDir','normal')
scatter3(obsx,obsy,obsz,1,'k')
set(gca, 'YAxisLocation', 'right')
text(min(xx)-dx(1), mean(yn),'$(b)$','interpreter', 'latex','FontSize',14);

% Add DO27 marker
scatter3(2000,3035,450,50,'*','r')
text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12,'BackgroundColor','w');

scatter3(2350,3600,450,50,'*','r')
text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12,'BackgroundColor','w');

%%
axes('Position',[0.58 .165 .38 .38])
temp =m2( padT+1:end-padB,padW+1:end-padE,padS+1:end-padN );
aa = isosurface(XX,YY,ZZ,temp,iso_cut);

% Create color map for depth
bb = interp1([max(zz);median(zz);min(zz)],[215 48 39;255 255 191;69 117 180]/255,aa.vertices(:,3));

view(cam_ang);
p = patch('Faces',aa.faces,'Vertices',aa.vertices,'FaceColor','b','LineStyle','none','FaceColor', 'flat','FaceVertexCData',bb);
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
hold on
% colorbar(gca,'NorthOutside')
axis equal
set(gca,'DataAspectRatio',[1 1 1/2])
camlight HEADLIGHT 
lighting phong
set(gca,'YDir','normal')
grid on
xlim([min(xx) max(xx)])
ylim([min(yy) max(yy)])
scatter3(obsx,obsy,obsz,1,'k')

% Add DO27 marker
scatter3(2000,3035,450,50,'*','r')

scatter3(2350,3600,450,50,'*','r')
%% Add colorbars
ax = axes('Position',[0.58 .215 .38 .38]);
bb = interp1([max(zz);median(zz);min(zz)],[215 48 39;255 255 191;69 117 180]/255,sort(aa.vertices(:,3)));
colormap(ax,bb);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',round([min(zz) median(zz) max(zz)]))
set(gca,'Visible','off');
text(0.35,1.25,'$Elevation (m)$', 'interpreter', 'latex','FontSize',14)

%% Plot Lp zones and alphas
lpmat = load([work_dir '\' norm_vec]);
[s,LP] = find_zones(lpmat);

mcell = nx*ny*nz;
% Smooth out the regions with 8-point averager
% Power determine the transition length
A = get_AVG_8pt(dx,dy,dz);

% A = A*A;
% A = spdiags(1./sum(A,2),0,mcell,mcell) *A;


t = A*(A*(A*s));

t = spdiags(1./sum(t,2),0,mcell,mcell) * t;

lp{1} = zeros(nx*ny*nz,1);
lp{2} = zeros(nx*ny*nz,1);
lp{3} = zeros(nx*ny*nz,1);
lp{4} = zeros(nx*ny*nz,1);

for jj = 1 : size(LP,1)
    
    lp{1}   =  lp{1} + t(:,jj)*LP(jj,1);
    lp{2}   =  lp{2} + t(:,jj)*LP(jj,2);
    lp{3}   =  lp{3} + t(:,jj)*LP(jj,3);
    lp{4}   =  lp{4} + t(:,jj)*LP(jj,4);

    
end

set(figure, 'Position', [50 0 775 1000]);

for jj = 1 : 4
m = reshape(lp{jj}(:),nz,nx,ny);

% Get top of model
m2D = get_model_top(m,nx,ny,nz,-100);

temp = reshape(m2D(padW+1:end-padE,padS+1:end-padN),nx-padE-padW,ny-padN-padS); 
temp = (temp');



if jj== 1
    ax = axes('Position',[0.1 (0.9 - 0.25) .35 .35]);
    h = imagesc(xx,yy,temp);
    caxis([0 2]);
    grid on
    text(min(xx)+dx(end)/2, min(yy)+dy(end)/2,'$l_p$','interpreter', 'latex','FontSize',14,'BackgroundColor','w');
    set(gca,'XTickLabel',[])
    ylabel('$y$', 'interpreter', 'latex','FontSize',14)
    set(gca,'YDir','normal')
    set(get(gca,'YLabel'),'Rotation',360);
    set(gca,'XTickLabel',[])
    set(gca,'YDir','normal')
elseif jj ==2
    ax = axes('Position',[0.5 (0.9 - 0.25) .35 .35]);
    h = imagesc(xx,yy,temp);
    caxis([0 2]);
    grid on
    text(min(xx)+dx(end)/2, min(yy)+dy(end)/2,'$l_{qx}$','interpreter', 'latex','FontSize',14,'BackgroundColor','w');
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])
    set(gca,'YDir','normal')
elseif jj==3
    ax = axes('Position',[0.1 (0.9 - 0.55) .35 .35]);
    h = imagesc(xx,yy,temp);
    caxis([0 2]);
    grid on
    text(min(xx)+dx(end)/2, min(yy)+dy(end)/2,'$l_{qy}$','interpreter', 'latex','FontSize',14,'BackgroundColor','w');
    ylabel('$y$', 'interpreter', 'latex','FontSize',14)
    set(gca,'YDir','normal')
    set(get(gca,'YLabel'),'Rotation',360);
    set(gca,'YDir','normal')
    xlabel('$x$', 'interpreter', 'latex','FontSize',14)
    
else 
    ax = axes('Position',[0.5 (0.9 - 0.55) .35 .35]);
    h = imagesc(xx,yy,temp);
    caxis([0 2]);
    grid on
    text(min(xx)+dx(end)/2,min(yy)+dy(end)/2,'$l_{qz}$','interpreter', 'latex','FontSize',14,'BackgroundColor','w');
    xlabel('$x$', 'interpreter', 'latex','FontSize',14)
    set(gca,'YTickLabel',[])
    set(gca,'YDir','normal')
    
end


axis square
hold on
caxis([0 2]);
hold on
scatter(obsx,obsy,1,'k')
% set(gca, 'YAxisLocation', 'right')


% Add DO27 marker
scatter3(2000,3035,450,50,'*','r')
text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12);

scatter3(2350,3600,450,50,'*','r')
text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12);



end

%% Add colorbars
ax = axes('Position',[0.325 0.05 .3 .3]);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 0.5 1])
set(cbar,'TickLabels',round([0 1 2]))
set(gca,'Visible','off');
text(0.4,1.35,'$norm$', 'interpreter', 'latex','FontSize',14)   

%% Plot Data
% Load data
[~, ~, ~, ~, ~, ~, ~, data, ~] = read_MAG3D_obs([work_dir '\' obsfile]);
% Load predicted l2
[~, ~, ~, ~, ~, ~, ~, dl2, ~] = read_MAG3D_obs([work_dir '\' predfile_l2]);
% Load predicted lp
[~, ~, ~, ~, ~, ~, ~, dlp, ~] = read_MAG3D_obs([work_dir '\' predfile_lp]);


nx = 100;
ny = 100;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / nx);
dy = (( ymax - ymin) / ny);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[X,Y] = meshgrid(x,y);

rl2 = (data - dl2)./wd;
rlp = (data - dlp)./wd;

data_interp     = griddata(obsx, obsy, data,X,Y,'linear'); 
% data_interp(isnan(data_interp)) = min(data_interp(:))*2;

dl2_interp     = griddata(obsx, obsy, dl2,X,Y,'linear'); 

dlp_interp     = griddata(obsx, obsy, dl2,X,Y,'linear');

rl2_interp        = griddata(obsx, obsy, rl2 ,X,Y,'linear');  

rlp_interp        = griddata(obsx, obsy, rlp ,X,Y,'linear'); 

set(figure, 'Position', [50 0 800 1000]);

axes('Position',[0.1 .68 .28 .28]);
h =imagesc(x,y,data_interp);hold on
set(h,'alphadata',~isnan(data_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
% set(gca,'XTickLabel',[])
grid on
axis equal
text(2100, 4600,'$Observed$', 'interpreter', 'latex','FontSize',14)
text(400, 2700,'$(a)$', 'interpreter', 'latex','FontSize',14)


axes('Position',[0.40 .68 .28 .28]);
h =imagesc(x,y,dl2_interp);hold on
set(h,'alphadata',~isnan(dl2_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
axis([min(x) max(x) min(y) max(y)])
set(gca,'YTickLabel',[])
grid on
axis equal
text(1600, 4600,'$Predicted\;(l_2-norm)$', 'interpreter', 'latex','FontSize',14)

axes('Position',[0.70 .68 .28 .28]);
h =imagesc(x,y,dlp_interp);hold on
set(h,'alphadata',~isnan(dlp_interp))
caxis([min(data_interp(:)) max(data_interp(:))]);
colormap(jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
set(gca,'YTickLabel',[])
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
axis([min(x) max(x) min(y) max(y)])
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
grid on
axis equal
text(1600, 4600,'$Predicted\;(l_p-norm)$', 'interpreter', 'latex','FontSize',14)

axes('Position',[0.40 .25 .28 .28]);
h =imagesc(x,y,rl2_interp);hold on
set(h,'alphadata',~isnan(rl2_interp))
caxis([-5 5])
colormap(jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
% xlabel('\bfEasting (m)')
ylabel('$y$', 'interpreter', 'latex','FontSize',14)
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(get(gca,'YLabel'),'Rotation',360);
axis([min(x) max(x) min(y) max(y)])
grid on
axis equal
text(3700, 4600,'$Residual$', 'interpreter', 'latex','FontSize',14)
text(400, 2700,'$(b)$', 'interpreter', 'latex','FontSize',14)

axes('Position',[0.70 .25 .28 .28]);
h =imagesc(x,y,rlp_interp);hold on
set(h,'alphadata',~isnan(rlp_interp))
caxis([-5 5])
colormap(jet);
scatter(obsx,obsy,2,'k.')
set(gca,'YDir','normal')
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
axis([min(x) max(x) min(y) max(y)])
xlabel('$x$', 'interpreter', 'latex','FontSize',14)
set(gca,'YTickLabel',[])
grid on
axis equal


%% Add colorbars
ax = axes('Position',[0.54 -.07 .28 .28]);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([-5 5]))
set(gca,'Visible','off');
text(0.3,1.35,'$\frac{(d^{pred} - d^{obs})}{\sigma} \;(nT)$', 'interpreter', 'latex','FontSize',12)

%% Add colorbars
ax = axes('Position',[0.54 0.36 .28 .28]);
cbar = colorbar('NorthOutside');
set(cbar,'Ticks',[0 1])
set(cbar,'TickLabels',round([min(data_interp(:)) max(data_interp(:))]))
set(gca,'Visible','off');
text(0.35,1.35,'$TMI (nT)$', 'interpreter', 'latex','FontSize',12)