function d = ampB(G,Wd,m)
% Function to compute amplitude of the field

d = zeros(size(G{1},1),1);
for ii = 1 : length(G)
    
    d = d + ( Wd * ( G{ii} * m ) ).^2 ;
    
end

d = sqrt(d);