function [s,p_vec] = find_zones(A)
% Function gets five vectors for the different regions of the model
% and output binary vectors

% A = [pvec qxvec qyvec qzvec rvec];
count = 1;
s = zeros(size(A,1),1);


% First find unique values
for ii = 1 : size(A,2)
    
    zi = unique(A(:,ii));
    
    for jj = 1 : length(zi)
        
        
        indx = A(:,ii)==zi(jj);
        
        % Check if the vector exist already
        tester = repmat(indx,1,size(s,2));
        check = tester == s;
        
        if sum(sum(check)==size(s,1)) ~= 0
            
            continue
        
        % Else, check that other vectors are constant
        else
            flag = 0;
            for kk = 1 : size(A,2)
              
                zone = A(indx,kk);
                
                if length(unique(zone))>1

                    flag = 1;

                end
            
            end
            
            if flag == 0
                
                s(:,count) = indx;
                count = count + 1;
                
            end
                
            
        end
        
    end
    
end

% Create lp vector for all zones
p_vec = zeros(1,5);
for ii = 1 : size(s,2)
    
    for jj = 1 : size(A,2)
        
       p_vec(ii,jj) = unique(A(s(:,ii)==1,jj));
       
    end
    
end
