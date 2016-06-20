% Create figure through 3D model
clear all
close all


addpath ..\..\FUNC_LIB;

%% Input Files
work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\MAG\Composite';

meshfile = 'Mesh_30m_padded.msh';

model2 = 'Composite_Tiled30m_l0l2.dat';
model1 = 'MAG3D_30m\maginv3d_012.sus';

% norm_vec = 'Lp_vec.txt';
% wmodel = 'Wvec.txt';

zpanel = 5;
ypanel = 10;

padE = 6;
padW = 6;

padN = 6;
padS = 6;

padT = 0;
padB = 12;

iso_cut = 0.2;
% obsfile = '..\..\DIGHEM_Mag_2pc_floor10nt_25m_ROT.obs';
% predfile_l2 = 'MAG3D_TMI_iter_2.pre';
% predfile_lp = 'MAG3D_TMI_iter_7.pre';

cam_ang = [30 60];

% Slices used for the animation. If only len=1 then only output a plot
% lvl = [1:10 9:-1:1];
lvl = [1:2:14];

maxcol = 0.1;
medcol = 0.05;

minc = 0.025;
intc = 0.025;
maxc = 0.2;


epsl = 1e-3;

ndv = -100;
%% Load in model and plot


% Load data
% [H, BI, BD, MI, MD, obsx, obsy, obsz, d, wd] = read_MAG3D_obs([work_dir '\' obsfile]);


[xn,yn,zn] = read_UBC_mesh([work_dir '\' meshfile]);

% Move to local coordinates
% obsx = obsx - xn(1);
% obsy = obsy - yn(1);
% obsz = obsz - zn(1);

%% Plot Data
% nx = 100;
% ny = 100;
% 
% xmin = ( min(obsx) );
% xmax = ( max(obsx) );
% ymin = ( min(obsy) );
% ymax = ( max(obsy) );
% 
% dx = (( xmax - xmin) / nx);
% dy = (( ymax - ymin) / ny);
% 
% x = xmin + cumsum(ones(1,nx)*dx);
% y = ymin + cumsum(ones(1,ny)*dy);
% 
% [X,Y] = meshgrid(x,y);
% 
% 
% data_interp     = griddata(obsx, obsy, d,X,Y,'linear'); 
% 
% % set(figure, 'Position', [50 0 800 1000]);
% figure;
% % axes('Position',[0.1 .68 .28 .28]);
% h =imagesc(x,y,data_interp);hold on
% set(h,'alphadata',~isnan(data_interp))
% caxis([min(data_interp(:)) max(data_interp(:))]);
% cmap = colormap(jet);
% 
% contour(x,y,data_interp,'k')
% set(gca,'YDir','normal')
% % xlabel('\bfEasting (m)')
% ylabel('$y$', 'interpreter', 'latex','FontSize',14)
% xlabel('$x$', 'interpreter', 'latex','FontSize',14)
% set(get(gca,'YLabel'),'Rotation',360);
% 
% quiver(1500,1500,250*(1/2),250*(sqrt(3)/2),'k','LineWidth',2,'MaxHeadSize',1);
% text(1500,1500,'$N$','interpreter', 'latex','FontSize',16,'Color','k','HorizontalAlignment','right','VerticalAlignment','top');
% 
% scatter(obsx,obsy,5,'k*')
% % Add DO27 marker
% scatter3(2000,3035,450,50,'o','r')
% text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12,'Color','w');
% 
% scatter3(2350,3600,450,50,'o','r')
% text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12,'Color','w');
%     
% % set(gca,'XTickLabel',[])
% grid on
% axis equal
% axis([min(x) max(x) min(y) max(y)])
% text(2100, 4600,'$Observed$', 'interpreter', 'latex','FontSize',14)
% text(400, 2700,'$(a)$', 'interpreter', 'latex','FontSize',14)

% Add color bar
    %%
%     ax = axes('Position',[0.68 0.25 .25 .5]);
%     colormap(ax,cmap);
%     cbar = colorbar('EastOutside');
%     set(cbar,'Ticks',[0 0.5 1])
%     set(cbar,'TickLabels',round([min(d) median(d) max(d)]*100)/100)
%     set(gca,'Visible','off');
%     text(1.2,1.1,'$(nT)$', 'interpreter', 'latex','FontSize',14)

%%
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

m2 = load([work_dir '\' model2]);
m2 = reshape(m2,nz,nx,ny);

set(figure, 'Position', [50 0 500 900]); 
for ii = 1: length(lvl);
% Get top of model
    m1_2D = get_model_top(m1,nx,ny,nz,ndv,lvl(ii));
    m1_2D = reshape(m1_2D(padW+1:end-padE,padS+1:end-padN),nx-padE-padW,ny-padN-padS); 
    m1_2D = (m1_2D');

    m2_2D = get_model_top(m2,nx,ny,nz,ndv,lvl(ii));

    m2_2D = reshape(m2_2D(padW+1:end-padE,padS+1:end-padN),nx-padE-padW,ny-padN-padS); 
    m2_2D = (m2_2D');
    
    ax = axes('Position',[0.2 .4 .75 .75]);
    h = imagesc(xx,yy,m1_2D);hold on
    set(h,'alphadata',m1_2D~=ndv)
    
    contour(xx,yy,m1_2D,minc:intc:maxc,'k')
    caxis([0 maxcol]);
    cmap =[1 1 1;0.900280117988586 0.963585436344147 0.990476191043854;0.800560235977173 0.927170872688293 0.980952382087708;0.700840353965759 0.89075630903244 0.971428573131561;0.601120471954346 0.854341745376587 0.961904764175415;0.501400589942932 0.817927181720734 0.952380955219269;0.401680678129196 0.78151261806488 0.942857146263123;0.301960796117783 0.745098054409027 0.933333337306976;0.328808456659317 0.754902005195618 0.897435903549194;0.355656117200851 0.764705896377563 0.861538469791412;0.382503777742386 0.774509787559509 0.82564103603363;0.40935143828392 0.7843137383461 0.789743602275848;0.436199098825455 0.79411768913269 0.753846168518066;0.463046759366989 0.803921580314636 0.717948734760284;0.489894419908524 0.813725471496582 0.682051301002502;0.51674211025238 0.823529422283173 0.64615386724472;0.543589770793915 0.833333373069763 0.610256433486938;0.570437431335449 0.843137264251709 0.574358999729156;0.597285091876984 0.852941155433655 0.538461565971375;0.624132752418518 0.862745106220245 0.502564132213593;0.650980412960052 0.872549057006836 0.466666668653488;0.677828073501587 0.882352948188782 0.430769234895706;0.704675734043121 0.892156839370728 0.394871801137924;0.731523394584656 0.901960790157318 0.358974367380142;0.75837105512619 0.911764740943909 0.32307693362236;0.785218715667725 0.921568632125854 0.287179499864578;0.812066376209259 0.9313725233078 0.251282066106796;0.838914036750793 0.941176474094391 0.215384617447853;0.865761697292328 0.950980424880981 0.179487183690071;0.892609357833862 0.960784316062927 0.143589749932289;0.919457018375397 0.970588207244873 0.107692308723927;0.946304678916931 0.980392158031464 0.0717948749661446;0.973152339458466 0.990196108818054 0.0358974374830723;1 1 0;1 0.966666638851166 0;1 0.933333337306976 0;1 0.899999976158142 0;1 0.866666674613953 0;1 0.833333313465118 0;1 0.800000011920929 0;1 0.766666650772095 0;1 0.733333349227905 0;1 0.699999988079071 0;1 0.666666686534882 0;1 0.633333325386047 0;1 0.600000023841858 0;1 0.566666662693024 0;1 0.533333361148834 0;1 0.5 0;1 0.466666668653488 0;1 0.433333337306976 0;1 0.400000005960464 0;1 0.366666674613953 0;1 0.333333343267441 0;1 0.300000011920929 0;1 0.266666680574417 0;1 0.233333334326744 0;1 0.200000002980232 0;1 0.16666667163372 0;1 0.133333340287209 0;1 0.100000001490116 0;1 0.0666666701436043 0;1 0.0333333350718021 0;1 0 0];
    colormap(cmap);
    axis square
    hold on
    grid on
    ylabel('$y$', 'interpreter', 'latex','FontSize',14)
%     xlabel('$x$', 'interpreter', 'latex','FontSize',14)
    set(gca,'YDir','normal')
    hold on
    set(gca,'XTickLabel',[])
    set(get(gca,'YLabel'),'Rotation',360);
%     scatter(obsx,obsy,1,'k')
    % set(gca, 'YAxisLocation', 'right')
    % text(min(xx)-dx(1), mean(yn),'$(a)$','interpreter', 'latex','FontSize',14);

    % Add DO27 marker
%     scatter3(2000,3035,450,50,'o','r')
    text(2000,2000,['$Depth: $' num2str(sum(dz(1:lvl(ii)))) ' m'],'interpreter', 'latex','FontSize',12);

%     scatter3(2350,3600,450,50,'o','r')
%     text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12);

%     text(1400,1400,'(a)', 'interpreter', 'latex','FontSize',14);
    %% Get top of model



    ax = axes('Position',[0.2 -.1 .75 .75]);
    h = imagesc(xx,yy,m2_2D); hold on
    set(h,'alphadata',m1_2D~=ndv)
    contour(xx,yy,m2_2D,minc:intc:maxc,'k')
    caxis([0 maxcol]);
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
    text(2000,2000,['$Depth: $' num2str(sum(dz(1:lvl(ii)))) ' m'],'interpreter', 'latex','FontSize',12);

%     scatter3(obsx,obsy,obsz,1,'k')
%     set(gca, 'YAxisLocation', 'right')
    % text(min(xx)-dx(1), mean(yn),'$(b)$','interpreter', 'latex','FontSize',14);

    % Add DO27 marker
%     scatter3(2000,3035,450,50,'o','r')
%     text(1600,3135,'$DO27$','interpreter', 'latex','FontSize',12);

%     scatter3(2350,3600,450,50,'o','r')
%     text(1950,3700,'$DO18$','interpreter', 'latex','FontSize',12);
%     text(1400,1400,'(b)', 'interpreter', 'latex','FontSize',14);
    %% Add color bar
    %%
    ax = axes('Position',[0.4 0.375 .4 .2]);
    colormap(ax,cmap);
    cbar = colorbar('NorthOutside');
    set(cbar,'Ticks',[0 0.5 1])
    set(cbar,'TickLabels',round([0 medcol maxcol]*100)/100)
    set(gca,'Visible','off');
    text(-.35,1.25,'$\kappa \;(SI)$', 'interpreter', 'latex','FontSize',14)

%     text(0.30,2,['Depth: ' num2str(sum(dz(3:lvl(ii))))],'BackgroundColor','w','FontSize',13,'EdgeColor','k')

    
    if length(lvl) > 1

        frame = getframe(figure(1));
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if ii == 1;
          imwrite(imind,cm,[work_dir '\field.gif'],'gif', 'Loopcount',inf,'DelayTime',1);
        else
          imwrite(imind,cm,[work_dir '\field.gif'],'gif','WriteMode','append','DelayTime',1);
        end

    end

end

%% Plot bar graph
% figure;
% axes('Position',[0.1 .2 0.75 0.6])
% [counts1, values1] = hist(m2(m2>-100 & m2 < 0.025), 150);
% [counts2, values2] = hist(m1(m1>-100 & m2 < 0.025), 150);
% bar(values1, counts1,'r','LineStyle','none');

% bar(values2, counts2,'b','LineStyle','none');

% xx= values1;
% Add epsilon value
% r = 1./(xx.^(2) + epsl.^2);
    
% dphi_dm = (xx).*r;

% [h_plot,h1,h2] = plotyy(values1,counts1,xx,dphi_dm,'bar','plot');  
% hold(h_plot(1),'on');
% hold(h_plot(2),'on');
% set(h2,'LineWidth',2,'Color','k')
% set(h1,'FaceColor',[0.65 0.65 0.65],'LineStyle','none')
% %         ylim(h_plot(1),[0 100]);
% xlim(h_plot(1),[0 0.025])
% xlim(h_plot(2),[0 0.025])
% set(h_plot(1),'yscale','log')
% set(h_plot(2),'ycolor','k')

% set(h_plot(2),'YTickLabel',[])
% axis(h_plot(1),'square')
% axis(h_plot(2),'square')

% bar(h_plot(1),values2, counts2,'FaceColor',[0.9 0.9 0.9],'LineStyle','none');
% plot(h_plot(2),[epsl epsl],[0 max(dphi_dm)],'k--')
% 
% h_plot(1).YTick = [1e+1 1e+4 1e+5];
% 
% % set(h_plot(2),'$\mathbf{\hat g_p}(m)$', 'interpreter', 'latex','FontSize',12)
% ylabel(h_plot(1),'$\mathbf{Hist}(m)$', 'interpreter', 'latex','FontSize',12)
% xlabel('m','interpreter', 'latex','FontSize',12)
% 
% ylabh = get(h_plot(1),'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') + [1e-3 0 0.00]);
%         
% text(0.0255,2.5,'$\mathbf{\hat g_p}(m)$','Rotation',-5,'LineWidth',1.5,'BackgroundColor','w','EdgeColor','k','interpreter', 'latex','FontSize',14,'HorizontalAlignment','left','VerticalAlignment','bottom')
% legend('S-IRLS','l2-norm')   
% text(epsl,0.8,'$\epsilon_p$','BackgroundColor','w','LineStyle','--','EdgeColor','k','interpreter', 'latex','FontSize',16,'HorizontalAlignment','left','VerticalAlignment','top')