function [m, tncg, Pac] = GN_PCG_solver( G, m, mref, nullcell, d, phi_in, beta , PreC, Pac, lowBvec, uppBvec , MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG )
% Solve H dm = J and update model

mcell = length(m);

Di = speye(mcell,mcell);

tncg = 0;        % Count total CG steps
solves = 1;       % Count number of GN solves
rddm = 1;       % Mesure the relative change in rel|dm|
ddm = [1 1];    % Mesure of change in model update |dm|

% sc_G = kron(spdiags([1;2;2],0,3,3),speye(mcell/3));

objfunc = @(x) sum( ( G * x - d ).^2 ) + ( x' * beta * MOF * x );
        A= [ G   ;...
            sqrt( beta ) * ( aVRWs ) ;...
            sqrt( beta ) * ( aVRWx ) ;...
            sqrt( beta ) * ( aVRWy ) ;...
            sqrt( beta ) * ( aVRWz )] ;
        
while solves < 5 && rddm > 1e-3

     switch FLAG

        case 'SMOOTH_MOD'

        b = [- (G * m - d) ; ...
            - sqrt( beta ) * ( aVRWs * (m - mref) ) ;...
            - sqrt( beta ) * ( aVRWx * m ) ;...
            - sqrt( beta ) * ( aVRWy * m ) ;...
            - sqrt( beta ) * ( aVRWz * m ) ];

        case 'SMOOTH_MOD_DIF'
            
        b = [- (G * m - d) ; ...
            - sqrt( beta ) * ( aVRWs * (m-mref) ) ;...
            - sqrt( beta ) * ( aVRWx * (m-mref) ) ;...
            - sqrt( beta ) * ( aVRWy * (m-mref) ) ;...
            - sqrt( beta ) * ( aVRWz * (m-mref) ) ];
     end
                    
                    

    %% Gauss-Newton step
    dm = zeros(mcell,1);
    [dm,~,ncg] = PCGLSQ( dm, A , b, PreC, Pac );
%     fprintf('CG iter %i \n',ncg);
    tncg = tncg + ncg;
    %% Step length, line search                

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
    phi_out = [phi_in phi_in];   
    while (phi_out(1) > phi_in || gamma == 2) && gamma > 1e-5 %phi_out(2)/phi_out(1) >= 1
        
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

    phi_in = phi_out(1);
    
    if solves == 1
        
        rddm = 1;
        ddm(1) = ddm(2);
        
    else
        
        rddm = ddm(2)/ddm(1);
    
    end


    % Update model
    m = m_temp;

    fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,rddm);
    
    solves = solves + 1;

end                