function [v] = Gvec(G,W,x)
% Function to compute the product of a cell-array of matrices with a vector

v = zeros(size(G{1},1),1);

if length(G) > 1
    for ii = 1 : length(G)

        v = v+W{ii} * (G{ii} * x);

    end
    
else
    
    v = W * (G{1} * x);
    
end