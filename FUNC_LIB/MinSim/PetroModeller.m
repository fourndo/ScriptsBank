% Function PetroModeller
%
% Program to calculate the statistics of a physical property model onto a
% geological block model

clear all
close all

%% INPUT PARAMETERS

addpath ..\
dsep = '\';


meshfile = '..\..\Mesh_20m.msh';
FMmodelfile = '..\..\Geology.dat';

XMLfile = '..\..\Geology.xml';

% % List of model files
% work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\Forward\Mesh_20m_Dens';
% % modelfile = 'Density_Randn_Std.den'; 
% modelfile = 'gzinv3d_001.den';
% ylbl{1} = 'Density Anomaly (g/cc)';
% modelOut = 'Density_Median_Std.den';
% axlim = ([-0.3 0.3]);
% ndv = -100000;


% % List of model files
% work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\Forward\Mesh_20m_Susc';
% % modelfile = 'maginv3d_001.sus'; ylbl{1} = 'Susceptibility';
% axlim = ([1e-5 1e-1]);
% ndv = -1;
% modelfile = 'maginv3d_002.sus'; ylbl{1} = 'Susceptibility';
% 
% work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\Forward\Mesh_20m_Cond';
% modelfile = 'Cond_Randn_Std_model.con'; ylbl{1} = 'Conductivity';
% modelOut = 'Median_Std_model.con';
% axlim = ([-5e-2 5e-2]);
axlim = ([-5 -0.1]);
% ndv = 1e-8;


work_dir = 'C:\Egnyte\Private\dominiquef\Projects\4414_Minsim\Modeling\Forward\Compilation\Models';
modelfile = 'Cond_Median_DH_Inv1D.con'; ylbl{1} = 'Conductivity';
modelOut = 'Cond_Median_Randn.con';
ndv = 1e-8;

dh_file = 'dhData';
% Outfile model
% modelOut = 'MedianSusc.sus';


[xn,yn,zn] = read_UBC_mesh([work_dir dsep meshfile]);

phymodel = load([work_dir dsep modelfile]);
FMmodel = load([work_dir dsep FMmodelfile]);

% load([work_dir dsep dh_file]);

% ndv = -100;
%% Load the classification
[ID,label] = Sort_XML([work_dir dsep XMLfile]);

% Replace underscores for spaces
ID = regexp(ID,'_','split');

% Define color scheme
ccode = [255 168 0;0 191 255;175 166 37;0 255 123;255 0 0;255 213 0;176 197 222;192 192 192]/255;
clabel = [{'Post-Intrusive'},{'Missi'},{'Undivided Intrusive'},{'Hidden'},{'VMS'},{'Millrock'},{'BlueLagoon'},{'Pre-Flin Flon'}];
cloc = [4;10;20;30;38;48;56;60];
%% Load both the physical property model and geology and compute stats
geoID = unique(FMmodel);

stat{1} = zeros(length(geoID),5);

for ii = 1 : length(geoID)
    
    
%     indx = dhData(:,5) == geoID(ii) & dhData(:,6) ~= 0 & dhData(:,6) ~= -99999;
    
%     if sum(indx)~=0
% %         % Density DH
% %         stat{1}(ii,1) = median(abs(dhData(indx,6))-2.7);
% %         stat{1}(ii,2) = std(abs(dhData(indx,6))-2.7);
% %         stat{1}(ii,3) = min(abs(dhData(indx,6))-2.7);
% %         stat{1}(ii,4) = max(abs(dhData(indx,6))-2.7);
% % 
% %         % Susc DH
%         stat{1}(ii,1) = median(abs(dhData(indx,6))/1000);
%         stat{1}(ii,2) = std(abs(dhData(indx,6))/1000);
%         stat{1}(ii,3) = min(abs(dhData(indx,6))/1000);
%         stat{1}(ii,4) = max(abs(dhData(indx,6))/1000);
%         
%     else
        
        indx = FMmodel == geoID(ii) & phymodel ~= ndv & phymodel ~= -9999;
%         
%         if sum(indx)~=0
%             
%             stat{1}(ii,1) = median(phymodel(indx));
%             stat{1}(ii,2) = std(phymodel(indx));
%             stat{1}(ii,3) = min(phymodel(indx));
%             stat{1}(ii,4) = max(phymodel(indx));
% 
%         end

        if sum(indx)~=0
            
            stat{1}(ii,1) = log10(median(phymodel(indx)));
            stat{1}(ii,2) = log10(std(phymodel(indx)));
            stat{1}(ii,3) = min(phymodel(indx));
            stat{1}(ii,4) = max(phymodel(indx));

        end
%     end
    
end

stat{1}(stat{1}(:,1)==0,1) = median(stat{1}(stat{1}(:,1)~=0,1));
stat{1}(stat{1}(:,2)==0,2) = median(stat{1}(stat{1}(:,2)~=0,2));

% load([work_dir dsep 'denStats.dat']);
% stat{1} = denStats(:,2:end);

% load([work_dir dsep 'magStats.dat']);
% stat{1} = magStats(:,2:end);
% 
%% Wisker Plot
set(figure, 'Position', [25 0 1200 1200]) 
zoff = 0;
lbl = 20;
for ii = 1
    

    subplot(1,3,ii)
        
    for jj = 1 : length(geoID);
        
        if stat{ii}(jj,1) == 0 || geoID(jj) == -9999
            
            continue
            
        end
        
        zoff = jj*30;
        yloc = (jj-2);

        cn = stat{ii}(jj,1);
%         ex = stat{ii}(jj,2);
        
        ex = abs(cn)*0.05;
        
        % Change color code for formation
        if geoID(jj) < 8
            
            rowc = ccode(1,:);
            
        elseif geoID(jj) >7 && geoID(jj) < 15
            
            rowc = ccode(2,:);
            
        elseif geoID(jj) >14 && geoID(jj) < 26
            
            rowc = ccode(3,:);
        elseif geoID(jj) >= 26 && geoID(jj) <= 36
            
            rowc = ccode(4,:);
            
        elseif geoID(jj) >= 37 && geoID(jj) <= 39
            
            rowc = ccode(5,:);
            
        elseif geoID(jj) >= 40 && geoID(jj) <= 54
            
            rowc = ccode(6,:);
            
        elseif geoID(jj) >= 55 && geoID(jj) <= 59
            
            rowc = ccode(7,:);
            
        else
            
            rowc = ccode(8,:);
            
        end
            
%         plot([stat{ii}(jj,3) stat{ii}(jj,4)],[yloc+zoff yloc+zoff],'k-','LineWidth',2); hold on
        h = fill([cn-ex cn-ex cn + ex cn + ex cn-ex],[yloc-lbl+zoff yloc+lbl+zoff yloc+lbl+zoff yloc-lbl+zoff yloc-lbl+zoff],rowc);
        plot(cn,yloc+zoff,'ko','MarkerSize',5,'MarkerFaceColor','k'); hold on
%         if kk == 1
%             temp = rname{rind(jj)};
%             text(cn,yloc+zoff,temp,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
%         end
% 
%         if sind(kk)==0
%             text(stat{ii}{jj}(kk,4),yloc+zoff,'Unclassified','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)
% 
%         else
        text(axlim(1),yloc+zoff,strcat(ID{label==geoID(jj)}{end-1:end}),'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10)
        
    end
    
   
    title(ylbl{ii})
    ylabel('Rock Unit');
    xlabel(ylbl{ii});

    colormap(jet)
%     ylim([0 25])
    grid on
    grid minor
end

% for ii = 1 : 8
%     
%     jj = cloc(ii);
%     zoff = jj*30;
%     yloc = (jj-2);
%     
%     h = fill([0.5 0.5 0.75 0.75 0.5],[yloc-lbl+zoff yloc+lbl+zoff yloc+lbl+zoff yloc-lbl+zoff yloc-lbl+zoff],ccode(ii,:));
%     text(0.5,yloc+zoff,clabel{ii},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)
% 
% end
set(gca,'YDir','reverse','YTickLabel',[])
subplot(1,3,1)
xlim(axlim)
% set(gca,'xscale','log');
ylim([0 1750])

%% Ouput a mean density model
m_out = ones(size(FMmodel))*ndv;

for ii = 2 : length(geoID)
    
    indx = FMmodel == geoID(ii);
    
%     rand_m = stat{1}(ii,1);

%     rand_m = randn(sum(indx),1)*stat{1}(ii,1)*0.05 + stat{1}(ii,1);
    rand_m = (randn(sum(indx),1))*10^(stat{1}(ii,1))*0.01 + 10^stat{1}(ii,1);

%     rand_m(rand_m<0) = 0;
    m_out(indx) = rand_m;
%     m_out(indx) = stat{1}(ii,1);
end
% 
save([work_dir dsep modelOut],'-ascii','m_out');