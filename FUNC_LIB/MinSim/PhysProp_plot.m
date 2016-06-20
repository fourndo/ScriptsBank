% Function Biplot_bar(X)
% Plotting script of geochem data for the CGS 
% The function creates a biplot with bar graphs OR pie charts
% for the elements of X.
%
% INPUT:
% X - Data matrix of size ndata-by-mvar
% It is assumed that the first two columns are the location on the biplot
% and all other columns are used to create the bar|pie graph.
%
% 
% Written by: D.Fournier
% Created : 15/07/2015

clear all
close all

%% INPUT PARAMETERS

work_dir = 'C:\LC\Private\dominiquef\Projects\4414_Minsim\Modeling\PhysProp\20160423';

% data_file = 'CHEM_den_Susc_VAR.txt';
litho_file = 'LithoCodes_Log.csv';
phys_file = 'LAS_combined_vDF.dat';
% phys_file{2} = 'log2.csv';
% phys_file{3} = 'log4.csv';
% data_file = 'CHEM_NRM_Susc_VAR.txt';

markfile = 'FM_markers.dat';

ndv = -99999;

% Add stratigraphy
strat_file = 'zone_geo.dat';

dsep = '\';

% Rock name
% rname = importdata([work_dir dsep 'Outfile.dat']);

% Column used for rock type and litho
rl = 1;
ll = 2;

% Flag for variable transformation: 0:Linear | 1:Log10
logflag = [0 0 0 0 0 1 0 0];

% Specify color scheme
cbar = [236 116 240;7 140 25;30 30 30;255 0 0;0 255 34]/255;

% Either bar | pie
var = 'pie';

% Specify bar thickness and length of bars as a fraction of the x-y axis
dx = 0.005;
dy = 0.075;

% OR if pie charts, radius of the pie
r=0.025;

ylbl{3} = 'Conductivity log(S/m)';
% axslim=[2.4 3.3 2.5 7.5];

ylbl{2} = 'Susceptibility log(SI)';
% xlbl = 'NRM log(A/m)';
% axslim=[-4 2 -4.5 0];

ylbl{1} = 'Density (g/cc)';

ccode = ['b', 'r', 'g' , 'm' , 'k', 'c']; 

lims = [2.4 3.6;0 300;-25 50;250 5000;3000 7500;1.5 5.5];
%% Load data and format
% Rcode is a from-to file
litho = importdata([work_dir dsep litho_file],',');
litho.litho = litho.textdata(2:end,:);

mark = importdata([work_dir dsep markfile]);
markname = mark.textdata(2:end,2);

%% Load database of physical properties   
dbase = importdata([work_dir dsep phys_file]);
dbase.label = dbase.textdata(2:end,1);

dhID = unique(lower(dbase.textdata(2:end,rl)));

%Re-sample by averaging segments
dz = 0.5;
dbase_dwns.data = [];
dbase_dwns.label = [];

for jj = 1 : length(dhID)
    
    indh = strcmpi(dhID{jj}, dbase.label(:,1));
    hphys = dbase.data(indh,:);
    
    endHole = max(hphys(:,4));
    data = [];
    label = [];
    count = 1;
    holeID{1} = dhID{jj};
    while dz*count < endHole
        
        % Find data within the interval
        indx = (hphys(:,4) >= (dz*count-dz/2)) & (hphys(:,4) < (dz*count+dz/2));
        
        if sum(indx)>0
            
            data = [data;mean(hphys(indx,:),1)];
            label = [label;holeID];            
        end
        
        count = count + 1;
        
    end
    
    dbase_dwns.data = [dbase_dwns.data;data];
    dbase_dwns.label = [dbase_dwns.label;label];
end
%% Find unique drill holes


count = 0;

undv{1} = 'ndv';
%% Cycle through the unique holes with phys prop and assign rock code
% Extract all data with same holeID
for jj = 1 : length(dhID)      

    % Get rock code hole
    indx = strcmpi(dhID{jj}, litho.litho(:)); 
    dh_intv = litho.data(indx,:);
    dh_lith = litho.litho(indx,:);
    
    % Get physprop interval for current hole
    indh = strcmpi(dhID{jj}, dbase_dwns.label(:,1));
    hphys = dbase_dwns.data(indh,:);
    
    rck_code = [];
    fm_code = [];
    
    for kk = 1 : sum(indh)% Cycle through the intervals

        indi = (hphys(kk,4) >= dh_intv(:,1)) & (hphys(kk,4) < dh_intv(:,2));

        % Make sure an interval is found
        if sum(indi) == 1
            
            rck_code = [rck_code;dh_lith(indi,5)];
            fm_code = [fm_code;dh_lith(indi,7)];
            
        elseif sum(indi) > 1
            
            ind = find(indi);
            rck_code = [rck_code;dh_lith(ind(1),5)];
            fm_code = [fm_code;dh_lith(ind(1),7)];
            
        else
            
            
            rck_code = [rck_code;undv];
            fm_code = [fm_code;undv];
            
        end

    end 

    % Append rock code to database
    dbase_dwns.label(indh,2) = rck_code(:);
    dbase_dwns.label(indh,3) = fm_code(:);
    
end


%% For each formation, plot all holes physical properties
ff = 7;
rcolm = 5;    
fcolm = 7;
intv = 50; % Interval around formation maker to plot

fmID = unique(lower(litho.litho(:,fcolm)));
cmap2 = hsv(length(fmID));

rockID = unique(lower(litho.litho(:,rcolm)));   
cmap1 = hsv(length(rockID));

%% Compute stats on rock type
for ii = 1 : size(dbase_dwns.data,2)-4
    
    % Get number of unique rock type
    rind = unique(dbase_dwns.label(:,2));
    
    for jj = 1 : length(rind);
        
        indx = strcmpi(rind(jj), dbase_dwns.label(:,2));
        
        hsub = dbase_dwns.data(indx,ii);
        lsub = dbase_dwns.label(indx,:);
        
        % Split again for formation
        sind = unique(lsub(:,3));        
        
        for kk = 1 : length(sind)
            
            ind =  strcmpi(lsub(:,3),sind(kk));
            
            stat{ii}{jj}(kk,1) = mean(hsub(ind),1);
            stat{ii}{jj}(kk,2) = std(hsub(ind),1);
            stat{ii}{jj}(kk,3) = min(hsub(ind),[],1);
            stat{ii}{jj}(kk,4) = max(hsub(ind),[],1);
            
        end
        
    end
    
end

%% Wisker Plot
set(figure, 'Position', [25 0 1200 1200]) 
zoff = 0;
lbl = 2;
for ii = 1 : size(dbase_dwns.data,2)-4
    

    subplot(1,3,ii)
    % Get number of unique rock type
    rind = unique(hphys{ii}(:,3));
    
    for jj = 2 : length(rind);
        
        indx = hphys{ii}(:,3)==rind(jj);
        
        hsub = hphys{ii}(indx,:);
        
        sind = unique(hsub(:,4));    
        
        for kk = 1 : length(sind);
            zoff = jj*30;
            yloc = (jj-2) + 4*kk ;
            
            cn = stat{ii}{jj}(kk,1);
            ex = stat{ii}{jj}(kk,2);

            plot([stat{ii}{jj}(kk,3) stat{ii}{jj}(kk,4)],[yloc+zoff yloc+zoff],'k-','LineWidth',2); hold on
            h = fill([cn - ex cn - ex cn + ex cn + ex cn - ex],[yloc-lbl+zoff yloc+lbl+zoff yloc+lbl+zoff yloc-lbl+zoff yloc-lbl+zoff],rind(jj)*2);

            if kk == 1
                temp = rname{rind(jj)};
                text(cn,yloc+zoff,temp,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
            end
            
            if sind(kk)==0
                text(stat{ii}{jj}(kk,4),yloc+zoff,'Unclassified','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)

            else
                text(stat{ii}{jj}(kk,4),yloc+zoff,stypeUnique{sind(kk)},'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',10)
            end
        end
        
    end
    
   
    title(ylbl{ii})
    ylabel('Rock Unit');
    xlabel(ylbl{ii});

    colormap(jet)
%     ylim([0 25])
    grid on
    grid minor
end
subplot(1,3,1)
xlim([2.5 3.5])

%% Cluster rock types based on physical properties


%% SCRIPT STARTS HERE

% % Load the data
% data = load([work_dir dsep data_file]);
% 
% ndata   = size(data,1);
% nvar    = size(data,2)-2;
% 
% % Data transfomation
% for ii = 1 : nvar
%     
%     if logflag(ii) == 1
%         
%         data(:,ii) = log10(data(:,ii));
%         
%     end
%     
%     
% end
% 
% % Create figure window
% set(figure, 'Position', [25 25 900 900]);
% 
% % Make biplot
% axes('Position',[0.075 .075 .9 .9]);
% scatter(data(:,1),data(:,2),5,'k'); hold on
% axis equal square
% axis(axslim)
% aratio = 1/1;
% 
% % Compute re-scaling parameter for graphing
% sclx = abs(axslim(2) - axslim(1));
% scly = abs(axslim(4) - axslim(3));
% 
% dx = dx * sclx;
% dy = dy * scly;
% % Random select some points
% % indx = randi(ndata,round(ndata/2),1);
% % 
% % data = data(indx,:);
% % ndata   = size(data,1);
% 
% 
% switch var
%     
%     case 'bar'
%         
%         % Thickness of bars
%         ddx = dx * abs( max(data(:,1)) - min(data(:,1)) );
% 
%         % Add bargraphs
%         for ii = 1 : ndata
% 
%             % Center of bar graph
%             xx = data(ii,1);
%             yy = data(ii,2);
% 
%             % Length of variables
%             lx = sum(data(ii,3:end));
% 
%             % Fraction of length
%             ddy = dy * data(ii,3:end) / lx;
% 
% 
%             for jj = 1 : nvar
% 
%                 y0 = yy - dy/2*lx + sum( ddy(1:(jj-1)));
% 
%                 patch([xx-ddx xx+ddx xx+ddx xx-ddx],...
%                     [y0 y0 y0+ddy(jj) y0+ddy(jj)],cbar(jj,:),'LineWidth',1.5);hold on
% 
% 
%             end
%             scatter(data(ii,1),data(ii,2),5,'k','filled');
% 
%         end
%     
%     case 'pie'
%         
%         % Add pie chart
%         for ii = 1 : ndata
% 
%             % Center of bar graph
%             xx = data(ii,1);
%             yy = data(ii,2);
% 
%             % Angles
%             ang = 0:0.1:pi;
%             nang = length(ang)*2;
% 
%             % Length of variables
%             rr =[r*sin(ang)* sclx;r*cos(ang) * scly];
% 
%             rr = [ [xx+rr(1,:);yy+rr(2,:)] [xx-rr(1,end:-1:1);yy+rr(2,end:-1:1)] ];
% 
% 
%             count = 1;
% 
%             for jj = 1 : nvar
% 
%                 % Fraction of length
%                 dpie = round( nang* data(ii,jj+2) / sum(data(ii,3:end)) );
% 
%                 if count + dpie > nang
% 
%                     dpie = nang- count;
% 
%                 end
% 
%                 patch([xx rr(1,count:count+dpie) xx],...
%                     [yy rr(2,count:count+dpie) yy],cbar(jj,:));hold on
% 
%                 count = count + dpie;
% 
%             end
%             scatter(data(ii,1),data(ii,2),5,'k','filled');
% 
%         end
% 
% end
% 
% 
% 
% xlabel(xlbl);
% ylabel(ylbl);