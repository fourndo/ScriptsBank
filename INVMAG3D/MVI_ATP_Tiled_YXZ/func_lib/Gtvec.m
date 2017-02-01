function [v] = Gtvec(G,W,x)
% Function to compute the product of a cell-array of matrices with a vector

v = zeros(size(G{1},2),1);

if length(G) > 1
    for ii = 1 : length(G)

        v = v+G{ii}' * (W{ii}' * x);

    end

else
    
    v = G{1}' * (W' * x);
    
end