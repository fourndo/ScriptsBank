% MagVector_Calc
% Function inputs a magnetization vector model and extract average
% declination and inclination of anomaly

clear all 
close all


addpath '..\.'
addpath 'C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\m_map';

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Paul_Lake\Data';

pole_file = 'Cande_1995_Normal_Polarity_TABLE.txt';
dsep = '\';

%% DATA TABLE
Pipes = ['Brent,Arnie,Mark,Hawk,Grizzly,Leslies,Zach,Roger,Jeager,Nora,Misery'];
Pipes = regexp(Pipes,',','split');
Ages = [47.1;47.5;47.5;48;50.8;52.1;52.8;67.6;69.1;2106;2480];
radio = [1 1 1 0 1 0 1 1 0];
invrt = [0 0 0 0 0 1 1 0 0 0 0];
%% Load polarity file
pole = load([work_dir dsep pole_file]);


set(figure, 'Position', [50 25 400 500]); 

ax1 = axes('Position',[0.05 .6 .9 .2]);
for ii = 1 : size(pole,1)
    
    patch([pole(ii,1) pole(ii,1) pole(ii,2) pole(ii,2)],[0 1 1 0],[0 0 0],'EdgeColor','none')
    hold on
    
end

axis([60 75 0 4])
set(gca,'Xdir','reverse','YTickLabel',[])
text(68,-1.5,'Time (Ma)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(73,-1,'Polarity','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(45,-1.5,'(b)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

for ii = 1 : size(Pipes,2)
    
    if mod(ii,3)==0
   
        quiver(Ages(ii,1),1.5,0,-0.5,0,'LineWidth',1.5,'Color','k')
   
    elseif ii ==1
        
        quiver(Ages(ii,1),3,0,-2.5,0,'LineWidth',1.5,'Color','k')
                
    elseif mod(ii,3)==1
           
        quiver(Ages(ii,1),2.5,0,-1.5,0,'LineWidth',1.5,'Color','k')
       
    else
           
   quiver(Ages(ii,1),2.0,0,-1.0,0,'LineWidth',1.5,'Color','k')
    end
    
end

for ii = 1 : size(Pipes,2)
    
    if mod(ii,3)==0
   text(Ages(ii,1),2, Pipes{ii},'HorizontalAlignment','center');
%    patch([1 1 0 0]/3+Ages(ii,1),[0 1 1 0]/3+1.5,[radio(ii) radio(ii) radio(ii)])
   
    elseif ii ==1
        text(Ages(ii,1),3.5, Pipes{ii},'HorizontalAlignment','center');
    elseif mod(ii,3)==1
           text(Ages(ii,1),3, Pipes{ii},'HorizontalAlignment','center');

    else
           text(Ages(ii,1),2.5, Pipes{ii},'HorizontalAlignment','center');
    end    
end

ax1 = axes('Position',[0.05 .48 0.025 0.025]);
patch([1 1 0 0],[0 1 1 0],[0 0 0])
text(0.5,2.5,'Polarity','interpreter', 'latex','FontSize',14)
text(1.1,0.2,'Normal','interpreter', 'latex','FontSize',14)
set(gca,'Visible','off');
axis equal

ax1 = axes('Position',[0.05 .44 0.025 0.025]);
patch([1 1 0 0],[0 1 1 0],[1 1 1])
text(1.1,0.2,'Reversed', 'interpreter', 'latex','FontSize',14)
set(gca,'Visible','off');
axis equal

%%
ax1 = axes('Position',[0.05 .2 .9 .2]);
for ii = 1 : size(pole,1)
    
    patch([pole(ii,1) pole(ii,1) pole(ii,2) pole(ii,2)],[0 1 1 0],[0 0 0],'EdgeColor','none')
    hold on
    
end

axis([45 55 0 4])
set(gca,'Xdir','reverse','YTickLabel',[])
% text(60,-1.5,'Time (Ma)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(73,3.5,'Polarity','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)
% text(45,-1.5,'(b)','HorizontalAlignment','center', 'interpreter', 'latex','FontSize',14)

for ii = 1 : size(Pipes,2)
    
    if mod(ii,3)==0
   
        quiver(Ages(ii,1),1.5,0,-0.5,0,'LineWidth',1.5,'Color','k')
   
    elseif ii ==1
        
        quiver(Ages(ii,1),3,0,-2.5,0,'LineWidth',1.5,'Color','k')
                
    elseif mod(ii,3)==1
           
        quiver(Ages(ii,1),2.5,0,-1.5,0,'LineWidth',1.5,'Color','k')
       
    else
           
   quiver(Ages(ii,1),2.0,0,-1.0,0,'LineWidth',1.5,'Color','k')
    end
    
end

for ii = 1 : size(Pipes,2)
    
    if mod(ii,3)==0
   text(Ages(ii,1),2, Pipes{ii},'HorizontalAlignment','center');
%    patch([1 1 0 0]/3+Ages(ii,1),[0 1 1 0]/3+1.5,[radio(ii) radio(ii) radio(ii)])
   
    elseif ii ==1
        text(Ages(ii,1),3.5, Pipes{ii},'HorizontalAlignment','center');
    elseif mod(ii,3)==1
           text(Ages(ii,1),3, Pipes{ii},'HorizontalAlignment','center');

    else
           text(Ages(ii,1),2.5, Pipes{ii},'HorizontalAlignment','center');
    end    
end
