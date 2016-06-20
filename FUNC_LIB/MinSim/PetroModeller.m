% Function PetroModeller
%
% Program to calculate the statistics of a physical property model onto a
% geological block model

clear all
close all

%% INPUT PARAMETERS

addpath ..\
dsep = '\';
work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\PhysProp\20160423';

meshfile = 'MinSim_VO_Local.msh';
FMmodelfile = 'Formation.dat';

XMLfile = 'Formation.xml';

% List of model files
modelfile = 'PYGRAV_lp.den'; ylbl{1} = 'Density';

% Outfile model
modelOut = 'MedianDensity.den';

[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);

phymodel = load([work_dir dsep modelfile]);
FMmodel = load([work_dir dsep FMmodelfile]);


ID = Sort_XML([work_dir dsep XMLfile]);
%% Load both the physical property model and geology and compute stats
geoID = unique(FMmodel);

stats{1} = zeros(length(geoID),5);

for ii = 1 : length(geoID)
    
    indx = FMmodel == geoID(ii) & phymodel ~= -100;
    
    if sum(indx)~=0

%         stats{1}(ii,3) = mean(phymodel(indx));
        stat{1}(ii,1) = median(phymodel(indx));
        stat{1}(ii,2) = std(phymodel(indx));
        stat{1}(ii,3) = min(phymodel(indx));
        stat{1}(ii,4) = max(phymodel(indx));
    end
    
end

%% Wisker Plot
set(figure, 'Position', [25 0 1200 1200]) 
zoff = 0;
lbl = 2;
for ii = 1
    

    subplot(1,3,ii)
        
    for jj = 1 : length(geoID);
        
        
        
        zoff = jj*30;
        yloc = (jj-2);

        cn = stat{ii}(jj,1);
        ex = stat{ii}(jj,2);

        plot([stat{ii}(jj,3) stat{ii}(jj,4)],[yloc+zoff yloc+zoff],'k-','LineWidth',2); hold on
        h = fill([cn - ex cn - ex cn + ex cn + ex cn - ex],[yloc-lbl+zoff yloc+lbl+zoff yloc+lbl+zoff yloc-lbl+zoff yloc-lbl+zoff],jj);

%         if kk == 1
%             temp = rname{rind(jj)};
%             text(cn,yloc+zoff,temp,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
%         end
% 
%         if sind(kk)==0
%             text(stat{ii}{jj}(kk,4),yloc+zoff,'Unclassified','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)
% 
%         else
        text(stat{ii}(jj,4),yloc+zoff,ID{geoID(jj)},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)
%         end
        
        
    end
    
   
    title(ylbl{ii})
    ylabel('Rock Unit');
    xlabel(ylbl{ii});

    colormap(jet)
%     ylim([0 25])
    grid on
    grid minor
end
% subplot(1,3,1)
% xlim([2.5 3.5])

%% Ouput a mean density model
m_out = zeros(size(phymodel));

for ii = 1 : length(geoID)
    
    indx = FMmodel == geoID(ii);
    m_out(indx) = stat{1}(ii,1);
    
end

save([work_dir dsep modelOut],'-ascii','m_out');