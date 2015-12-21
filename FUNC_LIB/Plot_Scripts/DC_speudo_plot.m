% Create dipole field from sphere and plot vectors

clear all 
close all

addpath ..\arrow3
addpath ..\gridfitdir

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\MtIsa\Modeling\DCIP2D';

datafile = 'obs_dc_PDP_East_edt.dat';
dsep = '\';

%% Load potential data

data = importdata([work_dir dsep datafile],' ',1);

data = data.data;


%% Plot the pseudo section
scz = cosd(45);

figure;
loc = data(:,1:4);
minx = min(loc(:));
maxx = max(loc(:));
minz = -(maxx - minx)/2 / scz ;
maxz = 0;



logd = sign(data(:,5)).*log10(abs(data(:,5)));



rx = mean(loc(:,3:4),2);
tx = loc(:,1);

midp = [mean([rx tx],2) -abs(rx - tx)/2 / scz];

set(figure(1), 'Position', [25 50 900 400])

for ii = 1 : size(data,1)
    axis equal
    axis([minx-50 maxx+50 minz maxz+50]); hold on
    colormap(jet)
    caxis([min(logd) max(logd)])
    axis([0 2000 -650 0])
    ylabel('Pseudo-Depth')
    xlabel('Easting')
%     title('Pseudo-section: Pole-Dipole')
    if loc(ii,1) == loc(ii,2)
       
        scatter(loc(ii,1),0,15,'rv'); hold on
        text(loc(ii,1),10,'Tx','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        scatter(loc(ii,3),0,15,'kv');
        scatter(loc(ii,4),0,15,'kv');
        text(mean(loc(ii,3:4)),10,'Rx','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','bottom')       
        
        scatter( midp(1:ii,1) , midp(1:ii,2) , 15 , logd(1:ii),'filled');
        plot([mean(loc(ii,3:4)) midp(ii,1)],[0 midp(ii,2)],'k:')
        plot([loc(ii,1) midp(ii,1)],[0 midp(ii,2)],'k:')
        
        hold off
    end
    set(gcf,'color','w');
    frame = getframe(figure(1));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ii == 1;
      imwrite(imind,cm,[work_dir '\field.gif'],'gif', 'Loopcount',inf,'DelayTime',0.05);
    else
      imwrite(imind,cm,[work_dir '\field.gif'],'gif','WriteMode','append','DelayTime',0.05);
    end

    clf
    
end

%% Grid data

xx = min(minx) : 10 : maxx;
zz = min(minz) : 10 : maxz;

[XX,ZZ] = ndgrid(xx,zz);

F = scatteredInterpolant(midp(:,1) , midp(:,2) , logd,'natural','none');

gridd = F(XX, ZZ);

for ii = 1 : 10

%     axes('Position',[0.05 0 .9 .9]);
    axis equal
    axis([minx-50 maxx+50 minz maxz+50]); hold on
    colormap(jet)
    caxis([min(logd) max(logd)])
    ylabel('Pseudo-Depth')
    xlabel('Easting')
    
    surface(XX,ZZ,gridd,'EdgeColor','none','FaceAlpha',ii/10);
    set(gcf,'color','w');    
    scatter( midp(:,1) , midp(:,2) , 5 , 'k','filled');
    axis([0 2000 -650 0])
    frame = getframe(figure(1));
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if ii == 10
        
        imwrite(imind,cm,[work_dir '\field.gif'],'gif','WriteMode','append','DelayTime',5);
        
    else
        
        imwrite(imind,cm,[work_dir '\field.gif'],'gif','WriteMode','append','DelayTime',0.1);
        
    end
    clf
end