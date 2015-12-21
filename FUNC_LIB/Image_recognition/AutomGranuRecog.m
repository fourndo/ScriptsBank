% Image recognition test

clear all 
close all

photo = imread('C:\Users\dominiquef.MIRAGEOSCIENCE\Dropbox\Master\FUNC_LIB\Image_recognition\gravel.jpg');
dphoto = double(photo);


% photo_bw = ind2grey(photo);
set(figure, 'Position', [10 250 1000 500]);
imagesc(photo)
axis equal

nX = size(photo,2);
nY = size(photo,1);

photo_bw = zeros(nY,nX);
 
dp_dx = (zeros(nY,nX));
dp_dy = (zeros(nY,nX));
dp_dxy = (zeros(nY,nX));
dp_dxy_avg = (zeros(nY,nX));

ddp_dx = (zeros(nY,nX));
ddp_dy = (zeros(nY,nX));

p_bw_avg = (zeros(nY,nX));
p_bw_var = (zeros(nY,nX));
p_bw_skw = (zeros(nY,nX));

gy =  (dphoto(:,:,1) + dphoto(:,:,2) + dphoto(:,:,3) ) /3;

% photo_bw(:,:,1) = gy;
% photo_bw(:,:,2) = gy;
photo_bw = gy;

set(figure, 'Position', [10 250 1000 500]);
imagesc(photo_bw);
axis equal

fft_photo_bw = fft2(photo_bw);

% F = fftshift(fft_photo_bw); % Center FFT
% 
% F = abs(F); % Get the magnitude
% F = log(F+1); % Use log, for perceptual scaling, and +1 since log(0) is undefined
% F = mat2gray(F); % Use mat2gray to scale the image between 0 and 1
% 
% imshow(F,[]); % Display the result

%% Window averaging
window=2;
for ii = (1+window) : (nX-window)
    
    for jj = (1+window) : (nY-window)
        

            p_bw_avg(jj,ii) = sum(sum(photo_bw((jj-window:jj+window),(ii-window:ii+window))))/9;
                    
%             p_bw_var(jj,ii) = (sum(sum((photo_bw((jj-window:jj+window),(ii-window:ii+window)) - p_bw_avg(jj,ii)).^2))/9);
%             
%             p_bw_skw(jj,ii) = sum(sum((photo_bw((jj-window:jj+window),(ii-window:ii+window)).^3)/p_bw_var(jj,ii)^3))/9;
%             
%             p_bw_kts(jj,ii) = sum(sum((photo_bw((jj-window:jj+window),(ii-window:ii+window)).^4)/p_bw_var(jj,ii)^4))/9;
    end

end
set(figure, 'Position', [10 250 1000 500]);
imagesc(p_bw_avg);title('window mean (dx)');
axis equal
hold on

trial = 5000;
hit=zeros(trial,2);

% Random selector for max
for ii=1:trial
    Xin = ceil(rand(1)*size(p_bw_avg,2));
    Yin = ceil(rand(1)*size(p_bw_avg,1));
    dp =1;

    if Xin>=size(p_bw_avg,2)-window

        Xin=Xin-window;

    elseif Yin>=size(p_bw_avg,1)-window

        Yin=Yin-window;
        
    elseif Xin<=window

        Xin=window+1;
        
    elseif Yin<=window

        Yin=window+1;

    end
    
while dp>0
    
    if Xin>window && Yin>window && Xin<=size(p_bw_avg,2)-window && Yin<=size(p_bw_avg,1)-window
        
        subset = p_bw_avg(Yin-window:Yin+window,Xin-window:Xin+window);
        
    elseif Xin>window && Yin<=window && Xin<=size(p_bw_avg,2)-window && Yin<=size(p_bw_avg,1)-window
        
        subset = p_bw_avg(Yin:Yin+window,Xin-window:Xin+window);
        
    elseif Xin<=window && Yin>window && Xin<=size(p_bw_avg,2)-window && Yin<=size(p_bw_avg,1)-window
        
        subset = p_bw_avg(Yin-window:Yin+window,Xin:Xin+window);
     
    elseif Xin>window && Yin>window && Xin>=size(p_bw_avg,2)-window && Yin<=size(p_bw_avg,1)-window

        subset = p_bw_avg(Yin-window:Yin+window,Xin-window:Xin);
        
    elseif Xin>window && Yin>window && Xin<=size(p_bw_avg,2)-window && Yin>=size(p_bw_avg,1)-window
        
        subset = p_bw_avg(Yin-window:Yin,Xin-window:Xin+window);
        
    end
    
    max_loc = max(max(subset));
    dp = max_loc - p_bw_avg(Yin,Xin);
    [Ystep,Xstep] = find(subset==max_loc);
    
    Xin = Xin + (Xstep(1)-(floor(size(subset,2)/2)+1));
    Yin = Yin + (Ystep(1)-(floor(size(subset,1)/2)+1));
  
end
    hit(ii,1:2) = [Xin Yin];
    plot(Xin,Yin,'r*')

end

figure(1);hold on; plot(hit(:,1),hit(:,2),'r*')
axis equal
% for ii = 1:10
%     
% %     group = p_bw_avg(hit(1,2)-ii:hit(1,2)+ii,hit(1,1)-ii:hit(1,1)+ii);
%     group = p_bw_avg(163-ii:163+ii,213-ii:213+ii);
%     
%     if ii ==1
%     group_mean = mean( mean(group) );
%     end
%     
%     N = (1+2*ii)^2;
%     
%     d = sqrt((group-group_mean).^2);
%     
%     CV(ii) = sqrt(1/N*sum(sum((group-group_mean).^2))) /  group_mean;
%     photo(163-ii:163+ii,213-ii:213+ii) = 1;
%     figure(1); imagesc(photo)
% end
%% First derivative
for ii = 1 : nX
    
    for jj = 1:nY
        

        if ii == nX 
            
           dp_dx(jj,ii) = p_bw_avg(jj,ii-1)-p_bw_avg(jj,ii);
           
        else
            
            dp_dx(jj,ii) = p_bw_avg(jj,ii+1)-p_bw_avg(jj,ii);
            
        end
        
            
    end

end
            

for ii = 1 : nX
    
    for jj = 1:nY
        

        if jj == nY 

           dp_dy(jj,ii) = p_bw_avg(jj-1,ii)-p_bw_avg(jj,ii);
           
        else

            dp_dy(jj,ii) = p_bw_avg(jj+1,ii)-p_bw_avg(jj,ii);
            
        end
        
            
    end

end



%%
for ii = 1 : nX
    
    for jj = 1:nY
        

        if jj == nY && ii == nX
            
            dp_dxy(jj,ii) = p_bw_avg(jj-1,ii-1,1) + p_bw_avg(jj,ii-1,1) + p_bw_avg(jj-1,ii,1) - 3*p_bw_avg(jj,ii,1);
          
        elseif jj == nY && ii == 1
            
            dp_dxy(jj,ii) = p_bw_avg(jj-1,ii+1,1) + p_bw_avg(jj,ii+1,1) + p_bw_avg(jj-1,ii,1) - 3*p_bw_avg(jj,ii,1);
          
        elseif jj == 1 && ii == 1
            
            dp_dxy(jj,ii) = p_bw_avg(jj+1,ii+1,1) + p_bw_avg(jj,ii+1,1) + p_bw_avg(jj+1,ii,1) - 3*p_bw_avg(jj,ii,1);
           
        elseif jj == 1 && ii == nX
            
            dp_dxy(jj,ii) = p_bw_avg(jj+1,ii-1,1) + p_bw_avg(jj,ii-1,1) + p_bw_avg(jj+1,ii,1) - 3*p_bw_avg(jj,ii,1);
         
        elseif  jj == nY && (ii ~= nX || ii ~= 1)
            
            dp_dxy(jj,ii) = p_bw_avg(jj-1,ii-1,1) + p_bw_avg(jj,ii-1,1) + p_bw_avg(jj-1,ii,1) + p_bw_avg(jj-1,ii+1,1) + p_bw_avg(jj,ii+1,1) - 5*p_bw_avg(jj,ii,1);
        elseif jj == 1 && (ii ~= nX || ii ~= 1)
            
            dp_dxy(jj,ii) = p_bw_avg(jj,ii-1,1) + p_bw_avg(jj+1,ii-1,1) + p_bw_avg(jj+1,ii,1) + p_bw_avg(jj,ii+1,1) + p_bw_avg(jj+1,ii+1,1) - 5*p_bw_avg(jj,ii,1);

        elseif ii == 1 && (jj ~= nY || jj ~= 1)
            
            dp_dxy(jj,ii) = p_bw_avg(jj+1,ii,1) + p_bw_avg(jj-1,ii,1) + p_bw_avg(jj-1,ii+1,1) + p_bw_avg(jj,ii+1,1) + p_bw_avg(jj+1,ii+1,1) - 5*p_bw_avg(jj,ii,1);

        elseif ii == nX && (jj ~= nY || jj ~= 1)
            
            dp_dxy(jj,ii) = p_bw_avg(jj-1,ii-1,1) + p_bw_avg(jj,ii-1,1) + p_bw_avg(jj+1,ii-1,1) + p_bw_avg(jj+1,ii,1) + p_bw_avg(jj-1,ii,1) - 5*p_bw_avg(jj,ii,1);

        else 
            
            dp_dxy(jj,ii) = p_bw_avg(jj-1,ii-1,1) + p_bw_avg(jj,ii-1,1) + p_bw_avg(jj+1,ii-1,1) + p_bw_avg(jj+1,ii,1) + p_bw_avg(jj-1,ii,1) + p_bw_avg(jj-1,ii+1,1) + p_bw_avg(jj,ii+1,1) + p_bw_avg(jj+1,ii+1,1) - 8*p_bw_avg(jj,ii,1);
          
            
        end
        
            
    end

end

%% Window averaging
for ii = 2 : nX-1
    
    for jj = 2:nY-1
        

            dp_dxy_avg(jj,ii) = (dp_dxy(jj,ii+1)+dp_dxy(jj+1,ii)+dp_dxy(jj,ii)+dp_dxy(jj,ii-1)+...
                dp_dxy(jj-1,ii)+dp_dxy(jj+1,ii+1)+dp_dxy(jj+1,ii-1)+dp_dxy(jj-1,ii+1)+dp_dxy(jj-1,ii-1))/9;
                    
            
    end

end

% dp_band1 = dp_dx_band1 + dp_dy_band1;

% figure; imagesc(dp_dx_band1);caxis([0 10])
% figure; imagesc(dp_dy_band1);caxis([0 10])
% figure; imagesc(dp_band1);caxis([0 10])

% dp_band2 = dp_dx_band2 + dp_dy_band2;

% figure; imagesc(dp_dx_band2);caxis([0 10])
% figure; imagesc(dp_dy_band2);caxis([0 10])
% figure; imagesc(dp_band2);caxis([0 10])

dp_tot_dxy = (dp_dx.^2 + dp_dy.^2).^0.5;

%% Second derivative
for ii = 1 : nX
    
    for jj = 1:nY
        

        if ii == nX 
            
           ddp_dx(jj,ii) = dp_tot_dxy(jj,ii-1)-dp_tot_dxy(jj,ii);
           
        else
            
            ddp_dx(jj,ii) = dp_tot_dxy(jj,ii+1)-dp_tot_dxy(jj,ii);
            
        end
        
            
    end

end
            

for ii = 1 : nX
    
    for jj = 1:nY
        

        if jj == nY 

           ddp_dy(jj,ii) = dp_tot_dxy(jj-1,ii)-dp_tot_dxy(jj,ii);
           
        else

            ddp_dy(jj,ii) = dp_tot_dxy(jj+1,ii)-dp_tot_dxy(jj,ii);
            
        end
        
            
    end

end


ddp_tot_dxy2 = (ddp_dx.^2 + ddp_dy.^2).^0.5;






% figure; imagesc(p_bw_var);title('window variance');
% 
% figure; imagesc(p_bw_skw);title('window skewness');caxis([0 1])
% figure; imagesc(p_bw_kts);title('window Kurtosis');caxis([0 1])

set(figure, 'Position', [10 250 1000 500]);
imagesc(dp_dx);title('horizontal derivative (dx)');%caxis([0 10])
axis equal

set(figure, 'Position', [10 250 1000 500]);
imagesc(dp_dy);title('vertical derivative (dy)');%caxis([0 10])
axis equal
% figure; imagesc(dp_dx_band3);caxis([0 10])
% figure; imagesc(dp_dy_band3);caxis([0 10])

set(figure, 'Position', [10 250 1000 500]);
imagesc(dp_tot_dxy);title('Total derivative) (dx+dy)');caxis([0 50])
axis equal




% figure; imagesc(dp_dxy_band1);
% figure; imagesc(dp_dxy_band2);
set(figure, 'Position', [10 250 1000 500]);
imagesc(ddp_tot_dxy2);title('Second order centered');caxis([-150 150])
axis equal

set(figure, 'Position', [10 250 1000 500]);
imagesc(dp_dxy_avg);title('First order centered AVG');caxis([-150 150])
axis equal