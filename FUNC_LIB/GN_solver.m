function [m, tncg, Pac] = GN_solver( G, J, m, mref, d, phi_in, beta , PreC, Pac, lowBvec, uppBvec , MOF, aVRWs)
% Solve H dm = J and update model

mcell = length(m);

Di = speye(mcell,mcell);
tol = 1e-5;      % Tolerance for the smallest step size possible
tncg = 0;        % Count total CG steps
solves = 1;       % Count number of GN solves
rddm = 1;       % Mesure the relative change in rel|dm|
ddm = [1 1];    % Mesure of change in model update |dm|

% sc_G = kron(spdiags([1;2;2],0,3,3),speye(mcell/3));

objfunc = @(x) sum( ( G(x(1),x(2)) - d ).^2 ) + ( x' * beta * MOF * x );
A = [ J; sqrt( beta ) * ( aVRWs )] ;
        
while solves < 5 && rddm > 1e-3

     

    b = [- (G(m(1),m(2)) - d) ; - sqrt( beta ) * ( aVRWs * (m - mref) ) ];

                    
                    

    %% Gauss-Newton step
    dm = zeros(mcell,1);
%     [dm,~,ncg] = PCGLSQ( dm, A , b, PreC, Pac );
    [dm,~,ncg] = CGiter(dm,A'*A,A'*b);
%     fprintf('CG iter %i \n',ncg);
    tncg = tncg + ncg;
    %% Step length, line search                

    temp = spdiags(Pac);
    
    % Combine active and inactive cells step if active bounds
    if sum(temp)~=mcell

        % Compute contribution from inactive cells
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

    % Reduce step length in order to reduce phi
    phi_out = [phi_in phi_in];   
    while (phi_out(1) > phi_in || gamma == 2) && gamma > tol %phi_out(2)/phi_out(1) >= 1
        
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

%     fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,rddm);
    
    solves = solves + 1;

end                