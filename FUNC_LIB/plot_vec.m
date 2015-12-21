function plot_vec(ah,X,Y,Z,U,V,W,val,tresh,dx,scale_m,scale_l,dwns)
% Plot vector with straight line and symbol scaled with length of vector


% dx = abs(X(1)-X(1,2))
if dwns~=1
    
    X = X(1:dwns:end,1:dwns:end);
    Y = Y(1:dwns:end,1:dwns:end);
    Z = Z(1:dwns:end,1:dwns:end);
    U = U(1:dwns:end,1:dwns:end);
    V = V(1:dwns:end,1:dwns:end);
    W = W(1:dwns:end,1:dwns:end);
    val = val(1:dwns:end,1:dwns:end);
    
end

X = X(:);
Y = Y(:);
Z = Z(:);
U = U(:);
V = V(:);
W = W(:);
val = val(:);

if isempty(Z)
    
    magvec = ( sqrt( U.^2 + V.^2 ) );
    vecsc = max(magvec);
    U = U .* (1 / vecsc +dx).^ scale_m ;
    V = V .* (1 / vecsc +dx).^ scale_m ;
else
    
    magvec = ( sqrt( U.^2 + V.^2 + W.^2 ) );
    U = U .* (magvec *dx).^ scale_m ;
    V = V .* (magvec *dx).^ scale_m ;
    W = W .* (magvec *dx).^ scale_m ;
end

if magvec ~= 0
    
%     U = sign(U).*(abs(U) / magvec * scale_m) .^ scale_l;
%     V = sign(V).*(abs(V) / magvec * scale_m) .^ scale_l;
%     W = sign(W).*(abs(W) / magvec * scale_m) .^ scale_l;

%     magvec = sqrt( U.^2 + V.^2 + W.^2 );

    
    
    
end

ii = val > tresh;

    
if isempty(Z)

    plot(ah,[X(ii)'; X(ii)'+U(ii)'],[Y(ii)'; Y(ii)'+V(ii)'],'k','LineWidth',1); hold on
    scatter(ah,X(ii),Y(ii),3,'k','filled');

else

    plot3(ah,[X(ii)'; X(ii)'+U(ii)'],[Y(ii)'; Y(ii)'+V(ii)'],[Z(ii)'; Z(ii)'+W(ii)'],'k','LineWidth',1); hold on
    scatter(ah,X(ii),Y(ii),3,'k','filled');
    %     scatter3(ah,X(ii),Y(ii),Z(ii),(val(ii)*sc + 1).^scale_l,val(ii),'filled');

   
end


        