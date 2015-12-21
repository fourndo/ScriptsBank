% Detrending data

close all
clear all

work_dir = 'C:\LC\Private\dominiquef\Projects\4239_Kaminak_Coffee_Mag\Modeling\Detrending\Part2';

data1 = load([work_dir '\14013_Coffee_NW_Infill_notie.dat']);
data2 = load([work_dir '\14013_Coffee_2_notie.dat']);
data3 = load([work_dir '\14013_Coffee_1_notie.dat']);
data4 = load([work_dir '\14013_Coffee_SE_Infill_notie.dat']);

reg1 = load([work_dir '\2011_Coffee_Regional_Mag_notie_Infill_out.dat']);
reg2 = load([work_dir '\2011_Coffee_Infill_1_Mag_notie.dat']);
reg3 = load([work_dir '\2011_Coffee_Infill_2_Mag_notie.dat']);

data = [data1;data2;data3;data4];
regional = [reg1;reg2;reg3]; regional(:,3) = regional(:,3)-57120;

dx=100;
dy=100;

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
X = reshape(X,nc,1);
Y = reshape(Y,nc,1);

F = scatteredInterpolant(data(:,1),data(:,2),data(:,3),'linear','none');
B = scatteredInterpolant(regional(:,1),regional(:,2),regional(:,3),'linear','none');
in_data = F(X,Y);
in_reg = B(X,Y);

keeper = isnan(in_data)==0;
keeper2 = isnan(in_reg)==0;
keeper = (keeper.*keeper2)==1;

in_data = in_data(keeper);
in_reg = in_reg(keeper);
X = X(keeper);
Y = Y(keeper);

cmin = min(in_data);
cmax = max(in_data);

residual = in_data - in_reg;
% Create system of linear equation
A = ones(length(residual),1);
A = [A X Y];
b = residual;

m = (A'*A)\(A'*b);

trend = (m(1) + X*m(2) + Y*m(3));
detrended = in_data - trend;


figure;scatter(X,Y,10,in_data);title('Local data');caxis([cmin cmax]);colorbar
figure;scatter(X,Y,10,in_reg);title('Regional data');caxis([cmin cmax]);colorbar
figure;scatter(X,Y,10,residual);title('Residual');colorbar
figure;scatter(X,Y,10,trend);title('Trend');colorbar
figure;scatter(X,Y,10,detrended);title('Original data minus trend');caxis([cmin cmax]);colorbar
