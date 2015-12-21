% function gridXY(x,y,data,x0,y0,dx,dy)
% Function takes an xydata and resample at location specified by the x0,y0
% and dx dy information

work_dir = 'C:\LC\Private\dominiquef\Projects\4239_Kaminak_Coffee_Mag\Modeling\GRIDDING';
data_file = '2011_Coffee_Infill_1_Mag_notie.dat';

dx = 25;
dy = 25;

data = load([work_dir '\' data_file]);

x0 = floor(min(data(:,1))); x0 = x0 - mod(x0,dx);
y0 = floor(min(data(:,2))); y0 = y0 - mod(y0,dy);

xmax = ceil(max(data(:,1)));
ymax = ceil(max(data(:,2)));

nx = ceil((xmax - x0)/dx);
ny = ceil((ymax - y0)/dy);
nc = nx * ny;

x = x0 + cumsum(ones(nx,1)*dx) - dx/2;
y = y0 + cumsum(ones(ny,1)*dy) - dy/2;

[X,Y] = ndgrid(x,y);

F = scatteredInterpolant(data(:,1),data(:,2),data(:,4));

interp_d = F(reshape(X,nc,1),reshape(Y,nc,1));

outfile = [reshape(X,nc,1) reshape(Y,nc,1) interp_d];

save([work_dir '\' data_file 'GRID.dat'],'-ascii','outfile');