% Function seek_max(obsx,obsy,data)
clear all
close all

work_dir = 'C:\Users\dominiquef.MIRAGEOSCIENCE\Google Drive\Research\Modelling\Amplitude\ES\I90D0';
obsfile = 'FWR_ES_lBl.dat';


[H, I, Dazm, D, obsx, obsy, obsz, data, ~] = read_MAG3D_obs([work_dir '\' obsfile]);

%% Create 2D matrix of interpolated data
nx = 200;
ny = 200;

xmin = ( min(obsx) );
xmax = ( max(obsx) );
ymin = ( min(obsy) );
ymax = ( max(obsy) );

dx = (( xmax - xmin) / nx);
dy = (( ymax - ymin) / ny);

x = xmin + cumsum(ones(1,nx)*dx);
y = ymin + cumsum(ones(1,ny)*dy);

[X,Y] = ndgrid(x,y);

% Grid data
data_interp     = griddata(obsx, obsy, data,X,Y,'natural'); 
data_interp     = flipud(data_interp);


set(figure, 'Position', [25 100 1800 900])
imagesc(x,y,data_interp);hold on
scatter(obsx,obsy,10,'*')
% xlabel('\bfEasting (m)')
ylabel('\bfNorthing (m)')
caxis([min(data) max(data)]);
axis equal


%% Seek peek anomaly using random seed generator
window=2;
trial = 100;
for ii=1:trial
    
    
    Xin = ceil(rand(1)*size(X,2));
    Yin = ceil(rand(1)*size(X,1));
    dp =1;

    if Xin>=size(X,2)-window

        Xin=Xin-window;

    elseif Yin>=size(X,1)-window

        Yin=Yin-window;
        
    elseif Xin<=window

        Xin=window+1;
        
    elseif Yin<=window

        Yin=window+1;

    end
    
while dp>0
    
    if Xin>window && Yin>window && Xin<=size(X,2)-window && Yin<=size(X,1)-window
        
        subset = X(Yin-window:Yin+window,Xin-window:Xin+window);
        
    elseif Xin>window && Yin<=window && Xin<=size(X,2)-window && Yin<=size(X,1)-window
        
        subset = X(Yin:Yin+window,Xin-window:Xin+window);
        
    elseif Xin<=window && Yin>window && Xin<=size(X,2)-window && Yin<=size(X,1)-window
        
        subset = X(Yin-window:Yin+window,Xin:Xin+window);
     
    elseif Xin>window && Yin>window && Xin>=size(X,2)-window && Yin<=size(X,1)-window

        subset = X(Yin-window:Yin+window,Xin-window:Xin);
        
    elseif Xin>window && Yin>window && Xin<=size(X,2)-window && Yin>=size(X,1)-window
        
        subset = X(Yin-window:Yin,Xin-window:Xin+window);
        
    end
    
    max_loc = max(max(subset));
    dp = max_loc - X(Yin,Xin);
    [Ystep,Xstep] = find(subset==max_loc);
    
    Xin = Xin + (Xstep(1)-(floor(size(subset,2)/2)+1));
    Yin = Yin + (Ystep(1)-(floor(size(subset,1)/2)+1));
  
end
    hit(ii,1:2) = [Xin Yin];
%     plot(Xin,Yin,'w*')
  
end

plot(hit(:,1),hit(:,2),'w*')
