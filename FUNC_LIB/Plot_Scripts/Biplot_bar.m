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

work_dir = '.\';

% data_file = 'CHEM_den_Susc_VAR.txt';
data_file = 'CHEM_den_POR_VAR.txt';
% data_file = 'CHEM_den_CHARG_VAR.txt';
% data_file = 'CHEM_den_RES_VAR.txt';
% data_file = 'CHEM_NRM_Susc_VAR.txt';

dsep = '\';

% Flag for variable transformation: 0:Linear | 1:Log10
logflag = [0 1 0 0 0 0 0 0];

% Specify color scheme
cbar = [236 116 240;7 140 25;30 30 30;255 0 0;0 255 34]/255;

% Either bar | pie
var = 'pie';

% Specify bar thickness and length of bars as a fraction of the x-y axis
dx = 0.005;
dy = 0.075;

% OR if pie charts, radius of the pie
r=0.025;

% Set axis limits
% ylbl = 'Susceptibility log(SI)';
% axslim=[2.4 3.3 -4.5 0];

ylbl = 'Porosity log(%)';
axslim=[2.4 3.3 -2.5 0.75];

% ylbl = 'Chargeability log(ms)';
% axslim=[2.4 3.3 -2.5 2.5];

% ylbl = 'Resistivity log(Ohms-m)';
% axslim=[2.4 3.3 2.5 7.5];

% ylbl = 'Susceptibility log(SI)';
% xlbl = 'NRM log(A/m)';
% axslim=[-4 2 -4.5 0];

xlbl = 'Density (g/cc)';

%% SCRIPT STARTS HERE

% Load the data
data = load([work_dir dsep data_file]);

ndata   = size(data,1);
nvar    = size(data,2)-2;

% Data transfomation
for ii = 1 : nvar
    
    if logflag(ii) == 1
        
        data(:,ii) = log10(data(:,ii));
        
    end
    
    
end

% Create figure window
set(figure, 'Position', [25 25 900 900]);

% Make biplot
axes('Position',[0.075 .075 .9 .9]);
scatter(data(:,1),data(:,2),5,'k'); hold on
axis equal square
axis(axslim)
aratio = 1/1;

% Compute re-scaling parameter for graphing
sclx = abs(axslim(2) - axslim(1));
scly = abs(axslim(4) - axslim(3));

dx = dx * sclx;
dy = dy * scly;
% Random select some points
% indx = randi(ndata,round(ndata/2),1);
% 
% data = data(indx,:);
% ndata   = size(data,1);


switch var
    
    case 'bar'
        
        % Thickness of bars
        ddx = dx * abs( max(data(:,1)) - min(data(:,1)) );

        % Add bargraphs
        for ii = 1 : ndata

            % Center of bar graph
            xx = data(ii,1);
            yy = data(ii,2);

            % Length of variables
            lx = sum(data(ii,3:end));

            % Fraction of length
            ddy = dy * data(ii,3:end) / lx;


            for jj = 1 : nvar

                y0 = yy - dy/2*lx + sum( ddy(1:(jj-1)));

                patch([xx-ddx xx+ddx xx+ddx xx-ddx],...
                    [y0 y0 y0+ddy(jj) y0+ddy(jj)],cbar(jj,:),'LineWidth',1.5);hold on


            end
            scatter(data(ii,1),data(ii,2),5,'k','filled');

        end
    
    case 'pie'
        
        % Add pie chart
        for ii = 1 : ndata

            % Center of bar graph
            xx = data(ii,1);
            yy = data(ii,2);

            % Angles
            ang = 0:0.1:pi;
            nang = length(ang)*2;

            % Length of variables
            rr =[r*sin(ang)* sclx;r*cos(ang) * scly];

            rr = [ [xx+rr(1,:);yy+rr(2,:)] [xx-rr(1,end:-1:1);yy+rr(2,end:-1:1)] ];


            count = 1;

            for jj = 1 : nvar

                % Fraction of length
                dpie = round( nang* data(ii,jj+2) / sum(data(ii,3:end)) );

                if count + dpie > nang

                    dpie = nang- count;

                end

                patch([xx rr(1,count:count+dpie) xx],...
                    [yy rr(2,count:count+dpie) yy],cbar(jj,:));hold on

                count = count + dpie;

            end
            scatter(data(ii,1),data(ii,2),5,'k','filled');

        end

end



xlabel(xlbl);
ylabel(ylbl);