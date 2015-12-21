% Image recognition test

clear all 
close all

photo = imread('C:\Users\dominiquef\Dropbox\DOM_Projects\Image_recognition\gravel.jpg');
dphoto = double(photo);


% photo_bw = ind2grey(photo);
figure; imagesc(photo)
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

figure;imagesc(photo_bw);

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
figure; imagesc(p_bw_avg);title('window mean (dx)');
hold on

trial = 5000;
hit=zeros(trial,2);

% Random selector for max
for ii=1:trial
    Xin = ceil(rand(1)*size(p_bw_avg,2));
    Yin = ceil(rand(1)*size(p_bw_avg,1));
    dp =1;

    if Xin==size(p_bw_avg,2)

        Xin=Xin-1;

    elseif Yin==size(p_bw_avg,1)

        Yin=Yin-1;

    end
    
while dp>0
    
    if Xin>1 && Yin>1 && Xin<size(p_bw_avg,2) && Yin<size(p_bw_avg,1)
        
        window = p_bw_avg(Yin-1:Yin+1,Xin-1:Xin+1);
        
    elseif Xin>1 && Yin==1 && Xin<size(p_bw_avg,2) && Yin<size(p_bw_avg,1)
        
        window = p_bw_avg(Yin:Yin+1,Xin-1:Xin+1);
        
    elseif Xin==1 && Yin>1 && Xin<size(p_bw_avg,2) && Yin<size(p_bw_avg,1)
        
        window = p_bw_avg(Yin-1:Yin+1,Xin:Xin+1);
     
    elseif Xin>1 && Yin>1 && Xin==size(p_bw_avg,2) && Yin<size(p_bw_avg,1)

        window = p_bw_avg(Yin-1:Yin+1,Xin-1:Xin);
        
    elseif Xin>1 && Yin>1 && Xin<size(p_bw_avg,2) && Yin==size(p_bw_avg,1)
        
        window = p_bw_avg(Yin-1:Yin,Xin-1:Xin+1);
        
    end
    
    max_loc = max(max(window));
    dp = max_loc - p_bw_avg(Yin,Xin);
    [Ystep,Xstep] = find(window==max_loc);
    
    Xin = Xin+Xstep(1)-2;
    Yin = Yin+Ystep(1)-2;
  
end
    hit(ii,1:2) = [Xin Yin];
    plot(Xin,Yin,'w*')
  
end

figure(1);hold on; plot(hit(:,1),hit(:,2),'w*')

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

figure; imagesc(dp_dx);title('horizontal derivative (dx)');%caxis([0 10])
figure; imagesc(dp_dy);title('vertical derivative (dy)');%caxis([0 10])
% figure; imagesc(dp_dx_band3);caxis([0 10])
% figure; imagesc(dp_dy_band3);caxis([0 10])
figure; imagesc(dp_tot_dxy);title('Total derivative) (dx+dy)');caxis([0 50])




% figure; imagesc(dp_dxy_band1);
% figure; imagesc(dp_dxy_band2);
figure; imagesc(ddp_tot_dxy2);title('Second order centered');caxis([-150 150])

figure; imagesc(dp_dxy_avg);title('First order centered AVG');caxis([-150 150])