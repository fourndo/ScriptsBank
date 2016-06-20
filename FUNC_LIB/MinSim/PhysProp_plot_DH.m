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
  

    % Cycle through the phys props and create new figure
    for jj = 1 : 6
        
        set(figure, 'Position', [25 0 2000 1200])
                
        count = 0;
        max_depth = -1;
        for ii = 1 : length(dhID)

            % Grab drillholes
            indx = strcmpi(dhID{ii}, dbase_dwns.label(:,1));
            
            sub_dbase.data = dbase_dwns.data(indx,:);
            sub_dbase.label = dbase_dwns.label(indx,:);
    
            % Check if Hidden Fm is intersected
            indx = find(strcmpi(fmID{ff}, sub_dbase.label(:,3)));

            if isempty(indx)==1
                
                continue
                
            else
                
                count = count + 1;
                
            end
            
            % Grab all data around the bottom of Hidden
            zmarker = sub_dbase.data(indx(end),4);
            indx = (sub_dbase.data(:,4) > zmarker-intv) & (sub_dbase.data(:,4) < zmarker+intv);
            
            % Create axis
            ax = axes('Position',[0.05+(count-1)*0.15 0.075 .1 .8]);
            sub_data = sub_dbase.data(:,jj+4);
            depth = sub_dbase.data(:,4);

            % Remove ndv
            indx2 = sub_data~=ndv;

            sub_data = sub_data(indx2);
            % Get the mean physprop for plotting
            x0 = mean([lims(jj,1) lims(jj,2)]);
            x = depth(indx2);
            
            % Shift depth to first marker
            x = x - zmarker;
            
            % Find the maximum depth for scaling the plot
%             if max_depth > min(x);
%                 
%                 max_depth = min(x);
%                 
%             end
            % Check if phys prop has to be logged
            if logflag(jj) == 1

                x = x(sub_data>0);
                sub_data = sub_data(sub_data>0);

                sub_data = log10(sub_data);

            end


            %plot(sub_data,x,'k'); hold on

    %         for kk = 1 : 2
    %             sub_data(2:end-1) = mean([sub_data(1:end-2) sub_data(3:end) sub_data(2:end-1)],2);
    %             
    %         end

            % Compute a running average
%             for kk = 6 : length(sub_data)-5
% 
%                 sub_data(kk) = mean(sub_data(kk-5:kk+5));
% 
%             end

            % Add rock markers
            indx = find(strcmpi(dhID{ii}, litho.litho(:,1)));
            
            for kk = 1 : length(indx)
               
                zz = litho.data(indx(kk),:) - zmarker;
                
                % Plot rock code
                plot([lims(jj,1) lims(jj,1)],[zz(1) zz(2)],'Color',...
                    cmap1(strcmpi(rockID,litho.litho(indx(kk),rcolm)),:),...
                    'LineWidth',10); hold on
                
                text(lims(jj,1),mean([zz(1) zz(2)]),litho.litho(indx(kk),rcolm),'HorizontalAlignment','center','VerticalAlignment','middle')
                
                % Plot formation code
                plot([lims(jj,2) lims(jj,2)],[zz(1) zz(2)],'Color',...
                    cmap2(strcmpi(fmID,litho.litho(indx(kk),fcolm)),:),...
                    'LineWidth',10); hold on
                
                text(lims(jj,2),mean([zz(1) zz(2)]),litho.litho(indx(kk),fcolm))
                
            end
            
            % Plot phys prop, mean and markers
            plot(sub_data,x,'r'); hold on
            plot(ones(length(x),1)*x0,x, 'k');

            indx3 = strcmpi(dhID{ii}, mark.textdata(2:end,1) );


%             plot(ones(sum(indx3),1)*x0, mark.data(indx3,4)- zmarker, 'ro','MarkerFaceColor','r');
%             text(ones(sum(indx3),1)*x0,mark.data(indx3,4)- zmarker,markname(indx3),'Rotation',45);
%             ax.YTickLabelRotation=90;
%             set(gca,'color','none')

            val = regexp(dbase.textdata(1,jj+5),'_','split');

            xlabel(val{1}(2));



            set(gca,'Ydir','reverse')

            if jj == 1

    %             set(gca,'YTickLabel',[], 'YTick',[]);
                %axis off
                ylabel('Depth (m)')
            end


    %         if jj ~= 1
    %             
    %             set(gca,'YTickLabel',[]);
    %             
    %         end

            if ~isempty(x)

                title(dhID{ii});

            end
            ylim([-100 100]);
            xlim([lims(jj,1) lims(jj,2)]);
            
        end

    end
% end
% %% Cycle through the unique holes with phys prop and assign rock code
% % Extract all data with same holeID
% for ii = 1 : size(dbase,2)
%     
%    %Get unique drill holes
%    dhUnique = unique(hphys{ii}(:,5));
%    
%    for jj = 1 : length(dhUnique)      
%        
%     % Get the size for pre-allocation
%     indx = strcmpi(dhID{dhUnique(jj)}, dh.textdata(:)); 
%     dh_in = dh.data(indx,:);
%     
%     indh = find(hphys{ii}(:,5) == dhUnique(jj));
%     
%         for kk = 1 : size(dh_in,1)% Cycle through the intervals
% 
%             indi = (dh_in(kk,1) < hphys{ii}(indh,1)) & (dh_in(kk,2) > hphys{ii}(indh,1)) & (hphys{ii}(indh,5) == dhUnique(jj));
% 
%             hphys{ii}(indh(indi),3) = dh_in(kk,3);
% 
%         end 
%         
%     end
%         
% %         dmat = [dmat;submat];
%         
% end
% 
% 
% %% REPEAT FOR FORMATION
% % Creating a matrix for all the rock intervals
% % Extract all data with same holeID
% 
% %% Load stratigraphy
% % Load stratigraphy file and reformat
% % if ~isempty(strat_file)
% %     
% % %     stype = load([work_dir dsep strat_file],'%s','%d','%d','%s');
% %     fid = fopen([work_dir dsep strat_file]);
% %     stype = textscan(fid, '%s%d%d%s','delimiter',',');
% %     fclose(fid);
% % end
% % 
% % % Get unique formation
% % stypeUnique = unique(lower(stype{4}(:)));
% % 
% % stype{5} = [stype{2}(:) stype{3}(:)];
% % for ii = 1 : length(stypeUnique)
% %     
% %     indx = strcmpi(stypeUnique{ii},stype{4}(:) );
% %     stype{5}(indx,3) = ii;
% %     
% % end
% % 
% % for ii = 1 : size(dbase,2)
% %     
% %    %Get unique drill holes
% %    dhUnique = unique(hphys{ii}(:,5));
% %    
% %    for jj = 1 : length(dhUnique)      
% %        
% %     % Get the size for pre-allocation
% %     indx = strcmpi(dhID{dhUnique(jj)}, stype{1}(:)); 
% %     dh_in = stype{5}(indx,:);
% %     
% %     indh = find(hphys{ii}(:,5) == dhUnique(jj));
% %     
% %         for kk = 1 : size(dh_in,1)% Cycle through the intervals
% % 
% %             indi = (dh_in(kk,1) < hphys{ii}(indh,1)) & (dh_in(kk,2) > hphys{ii}(indh,1)) & (hphys{ii}(indh,5) == dhUnique(jj));
% % 
% % %             if sum(indi)~=0 && ii == 1
% % %                 fprintf('stop')
% % %                 
% % %             end
% %             hphys{ii}(indh(indi),4) = dh_in(kk,3);
% % 
% %         end 
% %         
% %     end
% %         
% % %         dmat = [dmat;submat];
% %         
% % end
% 
% %% REPEAT FOR SECOND STRATIGRAPHY FOR UNCLASSIFIED SAMPLES
% % Add stratigraphy from SKUA model
% % Different format since only the downhole location of the interface is
% % provided. Need to convert to from-to interval
% 
% strat_file = 'SKUA_geology_mrkrs_try2.csv';
% 
% % Load stratigraphy file and reformat
% if ~isempty(strat_file)
%     
% %     stype = load([work_dir dsep strat_file],'%s','%d','%d','%s');
%     fid = fopen([work_dir dsep strat_file]);
%     stype = textscan(fid, '%s%d%s','delimiter',',');
%     fclose(fid);
% end
% 
% % Get unique formation
% stypeUnique2 = unique(lower(stype{3}(:)));
% 
% % Create from-to intervals here
% % Start with zero from and use the marker for lower limit
% stype{5} = [zeros(size(stype{2},1),1) stype{2}(:)];
%     
% %Get unique drill holes
% dhUnique = unique(stype{1}(:));
% 
% % Cycle drillholes and find multiple intervals, otherwise skip
% for ii = 1 : length(dhUnique)
%     
%     indx = find(strcmpi(dhUnique{ii},stype{1}(:) ));
%     
%     % If more than one interval on same hole, 
%     % get the interval from the one ahead
%     if length(indx) > 1
%         
%         stype{5}(indx(2:end),1) = stype{5}(indx(1:end-1),2);
%         
%     end
% 
%     
% 
% end
% 
% % Assign stratigraphy number to new 
% for ii = 1 : length(stypeUnique2)
%     
%     indx = strcmpi(stypeUnique2{ii},stype{3}(:) );
%     stype{5}(indx,3) = ii + size(stypeUnique,1);
%     
% end
% 
% stypeUnique = [stypeUnique;stypeUnique2];
% 
% for ii = 1 : size(dbase,2)
%     
%    %Get unique drill holes
%    dhUnique = unique(hphys{ii}(:,5));
%    
%    for jj = 1 : length(dhUnique)      
%        
%     % Get the size for pre-allocation
%     indx = strcmpi(dhID{dhUnique(jj)}, stype{1}(:)); 
%     dh_in = stype{5}(indx,:);
%     
%     indh = find(hphys{ii}(:,5) == dhUnique(jj));
%     
%         for kk = 1 : size(dh_in,1)% Cycle through the intervals
% 
%             indi = (dh_in(kk,1) < hphys{ii}(indh,1)) & (dh_in(kk,2) > hphys{ii}(indh,1)) & (hphys{ii}(indh,5) == dhUnique(jj));
% 
% %             if sum(indi)~=0 && ii == 1
% %                 fprintf('stop')
% %                 
% %             end
% 
%             % Only modify if unclassified:0
%             indu = hphys{ii}(indh,4) == 0;
%             
% %             if sum(indi & indu)~=0
% %                 sprintf('Stop!\n')
% %             end
%             
%             hphys{ii}(indh(indi & indu),4) = dh_in(kk,3);
% 
%         end 
%         
%     end
%         
% %         dmat = [dmat;submat];
%         
% end
% 
% 
% %% Compute stats on formations
% for ii = 1 : size(hphys,2)
%     
%     % Get number of unique formations type
%     fm_ind = unique(hphys{ii}(:,4));
%     
%     for jj = 1 : length(fm_ind);
%         
%         indx = hphys{ii}(:,4)==fm_ind(jj);
%         
%         hsub = hphys{ii}(indx,:);
%         
%         rock_ii = unique(hsub(:,3));        
%         
%         for kk = 1 : length(rock_ii)
%             
%             ind = hsub(:,3)==rock_ii(kk);
%             
%             stat{ii}{jj}(kk,1) = mean(hsub(ind,2),1);
%             stat{ii}{jj}(kk,2) = std(hsub(ind,2),1);
%             stat{ii}{jj}(kk,3) = min(hsub(ind,2),[],1);
%             stat{ii}{jj}(kk,4) = max(hsub(ind,2),[],1);
%             
%         end
%         
%     end
%     
% end
% 
% %% Wisker Plot
% set(figure, 'Position', [25 0 1500 1200]) 
% 
% 
% %% Start by plotting strat column
% axs{1} = axes('Position',[-.025 .1 .4 .8]);
% imag = imread([work_dir dsep 'Gibson_Strat_Column.PNG']);
% image(imag)
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% axis equal
% title('Stratigraphic Column - from Gibson et al. 2014')
% 
% %% Create sub plots for each formations and phys props and add wiskers
% zoff = 0;
% lbl = 0.5;
% for ii = 1 : size(hphys,2)
%     
% 
%     % Get number of unique rock type
%     fm_ind = unique(hphys{ii}(:,4));
%     
%     for jj = 2 : length(fm_ind);
%         
%         indx = hphys{ii}(:,4)==fm_ind(jj);
%         
%         hsub = hphys{ii}(indx,:);
%         
%         rock_ii = unique(hsub(:,3));    
%         
%         if fm_ind(jj) == 14 % Missi Fm
%            
%             axes('Position',[0.375+(ii-1)*0.215 .83 .2 .05]);
%             p_flag = 1;
%                         
% %         elseif fm_ind(jj) == 12 % Millrock Fm
% %            
% %             axes('Position',[0.375+(ii-1)*0.215 .32 .2 .18]);
% %             p_flag = 1;
% %             
% %         elseif fm_ind(jj) == 9 % Hidden Fm
% %            
% %             axes('Position',[0.375+(ii-1)*0.215 .51 .2 .22]);
% %             p_flag = 1;
% %             
% %         elseif  fm_ind(jj) == 4 % Blue Lagoon Fm
% %            
% %             axes('Position',[0.375+(ii-1)*0.215 .18 .2 .125]);
% %             p_flag = 1;
%             
%         elseif fm_ind(jj) == 20 % FlinFlon SKUA
%            
%             axes('Position',[0.375+(ii-1)*0.215 .1 .2 .41]);
%             p_flag = 1;
%             
%         elseif fm_ind(jj) == 21 % Hidden SKUA
%             
%             axes('Position',[0.375+(ii-1)*0.215 .51 .2 .22]);
%             p_flag = 1;
%             
%         elseif fm_ind(jj) == 22 % Louis SKUA
%             
%             axes('Position',[0.375+(ii-1)*0.215 .74 .2 .08]);
%             p_flag = 1;
%             
%         else
%             
%             p_flag = 0;
%             
%         end
%             
%         if p_flag==1    
%             for kk = 1 : length(rock_ii);
% 
%                 zoff = 0;
%                 yloc = kk ;
% 
%                 cn = stat{ii}{jj}(kk,1);
%                 ex = stat{ii}{jj}(kk,2);
% 
%                 plot([stat{ii}{jj}(kk,3) stat{ii}{jj}(kk,4)],[yloc+zoff yloc+zoff],'k-','LineWidth',2); hold on
%                 h = fill([cn - ex cn - ex cn + ex cn + ex cn - ex],[yloc-lbl+zoff yloc+lbl+zoff yloc+lbl+zoff yloc-lbl+zoff yloc-lbl+zoff],stat{ii}{jj}(kk,2));
% 
% %                 if kk == 1
% %                     temp = rname{fm_ind(jj)};
% %                     text(cn,yloc+zoff,temp,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',10)
% %                 end
% 
%                 if rock_ii(kk)==0
%                     text(stat{ii}{jj}(kk,1),yloc+zoff,'Unclassified','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10)
%                 else
%                     text(stat{ii}{jj}(kk,1),yloc+zoff,rname{rock_ii(kk)},'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10)
%                 end
%                 
%             end
%             
%             set(gca,'YTickLabel',[]);
%             
%             if fm_ind(jj) ~= 4 
%             
%                 set(gca,'XTickLabel',[]);
%                 
%             else
%                 
%                 xlabel(ylbl{ii})
%                 
%             end
%             if ii == 1
%                 xlim([2.6 3.6]);
%                 
%             elseif ii == 2
%                 
%                 xlim([-1 1.5]);
%                 
%             end
%             
% %             ylim = ([1.5 kk+0.5]);
% %             axis equal
%             caxis([0 0.25])
%             grid on
%             grid minor
%     
%         end
%     end
%     
%    
% 
% end
