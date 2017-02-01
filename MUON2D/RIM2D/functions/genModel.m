function [model] = genModel(nX, nZ, background, anomaly, target)
%Create density model using Xmax, Zmax,
%density of country rock, density of ore body
%Vector target gives the center position of the ore body
%Size of the target has been fixed to 2 meters across
model=ones(nZ, nX) .* background;


for ii = (target(1) + [-1 : 1])
    for kk = (target(2) + [-1 : 1])
    model(kk, ii) = anomaly;
    end
   
end

end
