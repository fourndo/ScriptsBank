function Laplacian_test(argin)
% Laplacian operator test
% [n , error] = Laplacian_test('type of mesh')
% This program computed the second order partial derivative of a function
% using the Laplacian operator. The test can be performed on a uniform
% staggered grid or on a tensor mesh. The script outputs a print-to-screen
% table of error between the analytical  and the computed solution.
%
% USER INPUT
% 'uniform': Perform the test on a uniform mesh
% 'tensor' : Perform the test on a padded mesh. The cell expansion factor 
% is computed every time the number cell is changed, so that the total
% dimensions of the grid is conserved.
% Function differentiated: 
% f(x) = cos(x/2).*sin(2*x)i + cos(x/2).*sin(2*x)j + cos(x/2).*sin(2*x)k
% 
% Analytical solution:
% V (f(x)) = -17/4*cos(x/2).*sin(2*x) - 2*sin(x/2).*cos(2*x)
%
%


    
    meshtype = argin;

% Perform the test for a sequence of n cells
for kk=2:5
    
    n=2^kk;   

    %% Create mesh for center and faces
    X0 = 0; Y0 = 0; Z0 = 0; %Origin
    
    switch meshtype 
    case 'uniform'
        if kk==2
            fprintf('**LAPLACIAN TEST**\nUniform mesh | | | | |\n')
            fprintf('Cell size\t|residual|\n')
        end  
        % %Uniform mesh
        % %Create mesh position (faces)
        xf = linspace(-pi,pi,n+1); xf=xf(:); 
        yf = linspace(-pi,pi,n+1); yf=yf(:);
        zf = linspace(-pi,pi,n+1); zf=zf(:);

        % Create mesh size
        dx = abs(xf(2:end) - xf(1:end-1));
        dy = abs(yf(2:end) - yf(1:end-1));
        dz = abs(zf(2:end) - zf(1:end-1));
    case 'tensor'
    %Tensor mesh
        if kk==2
            fprintf('**LAPLACIAN TEST**\nTensor mesh |    |   |  | ||\n')
            fprintf('Cell size\t|residual|\n')
        end  
        expanx = get_factor(n,0.1,pi);
        dx = [0.1*expanx .^[(n-1):-1:0]';0.1*expanx .^[0:(n-1)]'] ;
        % dx = [0.05*expanx .^[0:(n-1)]'] ;

        expany = get_factor((n),0.1,pi);
        dy = [0.1*expany .^[(n-1):-1:0]';0.1*expany .^[0:(n-1)]'] ;
        % dy = [0.05*expany .^[0:(n)]'] ;

        expanz = get_factor((n),0.1,pi);
        dz = [0.1*expanz .^[(n-1):-1:0]';0.1*expanz .^[0:(n-1)]'] ;
        % dz = [0.05*expanz .^[0:(n+1)]'] ;

        xf = [-pi;-pi+cumsum(dx)]; xf=xf(:); 
        yf = [-pi;-pi+cumsum(dy)]; yf=yf(:);
        zf = [-pi;-pi+cumsum(dz)]; zf=zf(:);
    otherwise
       
       fprintf('User forgot to specify a mesh type, or argument not recognized')
       fprintf('Two options: UNIFORM or TENSOR')
       fprintf('Laplacian test performed on UNIFORM mesh by default')
        xf = linspace(-pi,pi,n+1); xf=xf(:); 
        yf = linspace(-pi,pi,n+1); yf=yf(:);
        zf = linspace(-pi,pi,n+1); zf=zf(:);

        % Create mesh size
        dx = abs(xf(2:end) - xf(1:end-1));
        dy = abs(yf(2:end) - yf(1:end-1));
        dz = abs(zf(2:end) - zf(1:end-1));
    end
    



    % Create center-center mesh size (hmid)
    dxm = dx(1:end-1)/2 + dx(2:end)/2; dxm=[dxm(1);dxm;dxm(end)];
    dym = dy(1:end-1)/2 + dy(2:end)/2; dym=[dym(1);dym;dym(end)];
    dzm = dz(1:end-1)/2 + dz(2:end)/2; dzm=[dzm(1);dzm;dzm(end)];
    
    % Create hmid vectors
    dxn = dx(1:end-1)/2 + dx(2:end)/2; %dxn=[dx(1);dxn;dx(end)];
    dyn = dy(1:end-1)/2 + dy(2:end)/2; %dyn=[dy(1);dyn;dy(end)];
    dzn = dz(1:end-1)/2 + dz(2:end)/2;

        %Create mesh location (center)
    xc = X0 + dx(1)/2 + [0;cumsum(dxn)]; xc=xc(:); 
    yc = Y0 + dy(1)/2 + [0;cumsum(dyn)]; yc=yc(:);
    zc = Z0 + dz(1)/2 + [0;cumsum(dyn)]; zc=zc(:);
    

    % Compute number of faces and cells in 3D
    nxm = length(dxm); nym = length(dym) ; nzm = length(dzm);

    nx = length(dx); ny = length(dy) ; nz = length(dz);

    nfx = (nx+1) * (ny) * (nz);
    nfy = (nx) * (ny+1) * (nz);
    nfz = (nx) * (ny) * (nz+1);

    nface = nfx + nfy + nfz;

    mcell = nx * ny *nz;

    % Create diagonal matrices for hmid dimensions
    dXm = spdiags(1./dxm,0,nxm,nxm);
    dYm = spdiags(1./dym,0,nym,nym);
    dZm = spdiags(1./dzm,0,nzm,nzm);

    % Create cell-center dimension matrix
    dX = spdiags(1./dx,0,nx,nx); dXxl=spdiags([1./dx(1);1./dx;1./dx(end)],0,nx+2,nx+2);
    dY = spdiags(1./dy,0,ny,ny); dYxl=spdiags([1./dy(1);1./dy;1./dy(end)],0,ny+2,ny+2);
    dZ = spdiags(1./dz,0,nz,nz); dZxl=spdiags([1./dz(1);1./dz;1./dz(end)],0,nz+2,nz+2);

    Xif = kron( kron( ones(nz,1) , ones(ny,1) ), xf );
    Yif = kron( kron( ones(nz,1) , yf ), ones(nx,1) );
    Zif = kron( kron( zf ,ones(ny,1) ), ones(nx,1) );

    % nhx = (nx) * (ny+1) * (nz+1);
    % nhy = (nx+1) * (ny) * (nz+1);
    % nhz = (nx+1) * (ny+1) * (nz);

    ddx = @(n) spdiags (ones (n+1,1)*[-1,1],[0,1],n+1,n+2);

    %% Model parameters
    uo = 4*pi*10^-7;

    s_p = ones(mcell,1) * 1e-2; % Background conductivity [S/m]
    s_s = ones(mcell,1) * 1e-2; % Anomaly conductivity [S/m]


    %% Laplacian operator
    % Second partial Derivative in 3 directions
    %Partial derivatives for x-component
    d_dx=  dXxl.^0.5 * ddx(nx+1); d_dx=d_dx(:,2:end-1); dd_dx= dXm*d_dx'*d_dx; dd_dx([1,end],:)=0;
    DDxx = kron( kron( speye(nz) , speye(ny) ), dd_dx );

    d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY*d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
    DDyx = kron( kron( speye(nz) , dd_dy ), speye(nx+1) );

    d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ*d_dz'*d_dz; dd_dz([1,end],:)=0;
    DDzx = kron( kron( dd_dz , speye(ny) ), speye(nx+1) );

    % %%Partial derivatives for y-component
    d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX*d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
    DDxy = kron( kron( speye(nz) , speye(ny+1) ), dd_dx );

    d_dy= dYxl.^0.5 * ddx(ny+1); d_dy=d_dy(:,2:end-1); dd_dy= dYm*d_dy'*d_dy; dd_dy([1,end],:)=0;
    DDyy = kron( kron( speye(nz) , dd_dy ), speye(nx) );

    d_dz= dZm.^0.5 * ddx(nz); d_dz=d_dz(:,2:end-1); dd_dz= dZ*d_dz'*d_dz; dd_dz([1,end],:)=0;
    DDzy = kron( kron( dd_dz , speye(ny+1) ), speye(nx) );


    %%Partial derivatives for z-component
    d_dx= dXm.^0.5 * ddx(nx); d_dx=d_dx(:,2:end-1); dd_dx= dX*d_dx'*d_dx; dd_dx([1,end],:)=0;% dd_dx([1,end])=dd_dx([1,end])/2;
    DDxz = kron( kron( speye(nz+1) , speye(ny) ), dd_dx );

    d_dy= dYm.^0.5 * ddx(ny); d_dy=d_dy(:,2:end-1); dd_dy= dY*d_dy'*d_dy; dd_dy([1,end],:)=0;% dd_dy([1,end])=dd_dy([1,end])/2;
    DDyz = kron( kron( speye(nz+1) , dd_dy ), speye(nx) );

    d_dz= dZxl.^0.5 * ddx(nz+1); d_dz=d_dz(:,2:end-1); dd_dz= dZm*d_dz'*d_dz; dd_dz([1,end],:)=0;
    DDzz = kron( kron( dd_dz , speye(ny) ), speye(nx) );

    Oyx = sparse ( nfx , nfy );
    Ozx = sparse ( nfx , nfz );

    Oxy = sparse ( nfy , nfx );
    Ozy = sparse ( nfy , nfz );

    Oxz = sparse ( nfz , nfx );
    Oyz = sparse ( nfz , nfy );

    L_p = [DDxx+DDyx+DDzx Oyx Ozx;
         Oxy DDxy+DDyy+DDzy Ozy;%
         Oxz Oyz DDxz+DDyz+DDzz];


    % Face Variables
    [Xif, Yif, Zif] = ndgrid(xf, yc, zc);
    [Xjf, Yjf, Zjf] = ndgrid(xc, yf, zc);
    [Xkf, Ykf, Zkf] = ndgrid(xc, yc, zf);

%     % Function for the field
%     Eif = @(X,Y,Z) -Z .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
%     Ejf = @(X,Y,Z) -X .* Z .* exp(-5*(X.^2 + Y.^2 + Z.^2));
%     Ekf = @(X,Y,Z) -X .* Y .* exp(-5*(X.^2 + Y.^2 + Z.^2));
    Eif = @(X,Y,Z) sin(X);
    Ejf = @(X,Y,Z) sin(Y);
    Ekf = @(X,Y,Z) sin(Z);    
    Aif = Eif(Xif, Yif, Zif);
    Ajf = Ejf(Xjf, Yjf, Zjf);
    Akf = Ekf(Xkf, Ykf, Zkf);

        Af = [Aif(:); Ajf(:); Akf(:)];

%         % Attempt at a BC matrix
%         bc = zeros(nxm,1); bc(1)= f_x(xf(1)-dx(1))/dx(1)/dxm(1); bc(end)= f_x(xf(end)+dx(end))/dx(end)/dxm(end);
%         axx = kron(kron(ones(nz,1),ones(ny,1)),bc);
% 
%         bc = zeros(nym,1); bc(1)= f_x(yf(1)-dy(1))/dy(1)/dym(1); bc(end)= f_x(yf(end)+dy(end))/dy(end)/dym(end);
%         ayy = kron(kron(ones(nz,1),bc),ones(nx,1));
% 
%         bc = zeros(nzm,1); bc(1)= f_x(zf(1)-dz(1))/dz(1)/dzm(1); bc(end)= f_x(zf(end)+dz(end))/dz(end)/dzm(end);
%         azz = kron(kron(bc,ones(ny,1)),ones(nx,1));
% 
%         BC = [axx;ayy;azz];
        LAf = L_p*Af;
%         LAfBC = LAf -BC;

    %     dd_fx = @(x) -sin(x) ;
%         dd_fx = @(x,y,z) (70 .* z .* y - 100 * ( z .* y .* x.^2 + z .* y.^3 + z.^3 .* y)) .* exp(-5*(x.^2 + y.^2 + z.^2));
%         dd_fy = @(x,y,z) (70 .* z .* x - 100 * ( z .* y .* y.^2 + z .* x.^3 + z.^3 .* x)) .* exp(-5*(x.^2 + y.^2 + z.^2));
%         dd_fz = @(x,y,z) (70 .* x .* y - 100 * ( x .* y .* z.^2 + x .* y.^3 + x.^3 .* y)) .* exp(-5*(x.^2 + y.^2 + z.^2));
        dd_fx = @(x,y,z) -sin(x);
        dd_fy = @(x,y,z) -sin(y);
        dd_fz = @(x,y,z) -sin(z);
        
        Bif = dd_fx(Xif,Yif,Zif);
        Bjf = dd_fy(Xjf,Yjf,Zjf);
        Bkf = dd_fz(Xkf,Ykf,Zkf);

    %     Bkf = -(15/4)*cos(Zif/2).*sin(2*Zif) - sin(Zif/2).*cos(2*Zif);

        Bf = [Bif(:); Bjf(:); Bkf(:)];

    %    fprintf('\t1/h\t\t|f(x)-f(x+h)|\t|f(x+h) - f(x)- dfdx*h|\n') 
       residual(kk,1)= n; 
       residual(kk,2)= norm(abs(LAf) - abs(Bf), 'inf');

       fprintf('%3.2e   %3.2e  \n',1/n,residual(kk,2))

end
   fprintf('**end of test**\n\n')