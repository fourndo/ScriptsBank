function [m,tncg,Pac] = GNsolver( G, m, d, phi, MOF, PreC, Pac, beta , aVRWs, aVRWx, aVRWy, aVRWz )
% Solve H dm = J and update model

mcell = length(m);

tncg = 0;        % Count total CG steps
solves = 1;       % Count number of GN solves
rddm = 1;       % Mesure the relative change in rel|dm|
ddm = [1 1];    % Mesure of change in model update |dm|
Di = speye(mcell,mcell);


lowBvec = ones(mcell,1) * -1;
uppBvec = ones(mcell,1) * 1;

objfunc = @(x) sum( ( G * x - d ).^2 ) + ( x' * beta * MOF * x );
cg_count = 0;
while solves < 10 && rddm > 1e-3

    if isempty(aVRWz)==0

        A= [ G  ;...
            sqrt( beta ) * ( aVRWs ) ;...
            sqrt( beta ) * ( aVRWx ) ;...
            sqrt( beta ) * ( aVRWy ) ;...
            sqrt( beta ) * ( aVRWz )] ;

        b = [- (G *m - d) ; ...
            -sqrt( beta ) * ( aVRWs * m) ;...
            -sqrt( beta ) * ( aVRWx * m);...
            -sqrt( beta ) * ( aVRWy * m);...
            -sqrt( beta ) * ( aVRWz * m)] ;


    elseif isempty(aVRWy)==0

        A= [ G  ;...
            sqrt( beta ) * ( aVRWs ) ;...
            sqrt( beta ) * ( aVRWx ) ;...
            sqrt( beta ) * ( aVRWy )] ;

        b = [- (G *m - d) ; ...
            -sqrt( beta ) * ( aVRWs * m) ;...
            -sqrt( beta ) * ( aVRWx * m);...
            -sqrt( beta ) * ( aVRWy * m)] ;

    else


        A= [ G  ;...
            sqrt( beta ) * ( aVRWs ) ;...
            sqrt( beta ) * ( aVRWx )] ;

        b = [- (G *m - d) ; ...
            -sqrt( beta ) * ( aVRWs * m) ;...
            -sqrt( beta ) * ( aVRWx * m)] ;


    end

    %% Gauss-Newton step
    dm = zeros(mcell,1);
%                 [dm,~,ncg] = CGLSQ( dm, A , b );
    [dm,~,ncg] = PCGLSQ( dm, A , b, PreC, Pac );
    tncg = tncg + ncg;
    %% Step length, line search                
%                 ncg = ncg+iter; % Record the number of CG iterations

    temp = spdiags(Pac);

     % Combine active and inactive cells step if active bounds
    if sum(temp)~=mcell

        rhs_a = ( Di - Pac ) * (A' * b);
        dm_i = max( abs( dm ) );
        dm_a = max( abs(rhs_a) ); 

        if dm_i < dm_a
            dm = dm + rhs_a * dm_i / dm_a /5 ;
        else
            dm = dm + rhs_a;
        end

    end

    gamma = 2;

    % Reduce step length in order to reduce phid
    phi_out = [phi phi];    
    while (phi_out(1) > phi || gamma == 2) && phi_out(2)/phi_out(1) >= 1

        phi_out(2) = phi_out(1);

        gamma = 0.5 * gamma;

        gdm = gamma * dm;

        ddm(2) = norm(gdm);

        m_temp = m + gdm;

        lowb = m_temp <= lowBvec;
        uppb = m_temp >= uppBvec;

        % Apply bound on model
        m_temp(lowb==1) = lowBvec(lowb==1);
        m_temp(uppb==1) = uppBvec(uppb==1);

        % Update projection matrix
        Pac = spdiags((lowb==0).*(uppb==0),0,mcell,mcell);

        phi_out(1) = objfunc(m_temp);

    end

    % Mesure relative change in dm
    if solves == 1

        rddm = 1;
        ddm(1) = ddm(2);

    else

        rddm = ddm(2)/ddm(1);

    end


    % Update model
    m = m_temp;

%     fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,rddm);
%     fprintf('CG iter %i \n',ncg);
    solves = solves + 1;

end                