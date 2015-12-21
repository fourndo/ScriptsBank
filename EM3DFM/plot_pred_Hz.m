function [d_pred] = plot_pred_Hz(work_dir,obsfile)

fid = fopen([work_dir '\' obsfile],'r');

% ndata = size(dwndata.data,1);
line = fgets(fid);
d_pred = zeros(1,16);
freq = [];
count = 1;
while line~=-1
    
    if isempty(regexp(line,'frequency','match'))==0
        
        temp = regexp(line,'[=]','split');
        freq = str2double(temp{2});
        
    end
        
    if isempty(regexp(line,'%','match'))==1 && isempty(regexp(line,'\d\s','match'))==0
        
        d_pred(count,1) = freq;
        d_pred(count,2:end) = str2num(line);
        count = count+1;
        
    end
    
    line = fgets(fid);
    
    
end

fclose(fid);

freq = unique(d_pred(:,1));

% Grab the XY of first frequency
X = d_pred(d_pred(:,1)==freq(1),2); 
Y = d_pred(d_pred(:,1)==freq(1),3);

% Plot data
xmin = floor( min(X) );
xmax = ceil( max(X) );
ymin = floor(  min(Y) );
ymax = ceil(  max(Y) );

dx = (xmax-xmin)/100;
dy = (ymax-ymin)/100;

nx = ceil(( xmax - xmin) / dx);
ny = ceil(( ymax - ymin) / dy);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[Xgrid,Ygrid] = meshgrid(x,y);

%% Plot Frequencies
% Obs_R_interp = griddata(X, Y, dwndata.getData(1:ndata/3,37),Xgrid,Ygrid); 
 figure;
for ii = 1 : length(freq);
    
    index = d_pred(:,1)==freq(ii);
    real_interp = griddata(X, Y, d_pred(index,15),Xgrid,Ygrid,'cubic');
    imag_interp = griddata(X, Y, d_pred(index,16),Xgrid,Ygrid,'cubic');

   
    subplot(3,2,2*(ii-1)+1)
    pcolor(x,y,real_interp); shading interp; 
    title(['Real - ' num2str(freq(ii)) 'Hz']); 
    h = colorbar;
    xlabel('Easting','fontsize',14); 
    ylabel('Northing','fontsize',14); 
    axis equal

    subplot(3,2,2*(ii-1)+2)
    pcolor(x,y,imag_interp); shading interp; 
    title(['Imag - ' num2str(freq(ii)) 'Hz']); 
    h = colorbar;
    xlabel('Easting','fontsize',14); 
    ylabel('Northing','fontsize',14); 
    axis equal

end