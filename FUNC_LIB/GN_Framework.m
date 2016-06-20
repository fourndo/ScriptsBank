% Framework for Gauss-Newton solves
% D. Fournier
% 13-01-16

% Save some inversion parameters
phi     = []; % Obj Func
phi_m   = []; % Model Obj Func
phi_d   = []; % Misfit Func
dphim   = []; % Change in Model Obj Func
beta    = []; % Trade-off parameter
Pac     = []; % Projection matrix for active cells


while dm > tol
    
%% Jacobi Pre-conditionner
    diagA = sum(G.^2,1) + beta(count)*spdiags(MOF,0)';
    PreC     = Pac * spdiags(1./diagA(:),0,mactv,mactv);

    %% Gauss-Newton steps

    fprintf('\n# # # # # # # # #\n');
    fprintf('BETA ITER: \t %i  \nbeta: \t %8.5e \n',count,beta(count));
    fprintf('eps_p: \t %8.5e \t eps_p*: \t %8.5e\n',delta_p(count),eps_p)
    fprintf('eps_q: \t %8.5e \t eps_q*: \t %8.5e\n',delta_q(count),eps_q)
    % Save previous model before GN step
    
    
    [invmod, ncg, Pac] = GN_PCG_solver( G, invmod, mref, nullcell, d, phi(end), beta(count) , PreC, Pac, lowBvec, uppBvec, MOF, aVRWs, aVRWx, aVRWy, aVRWz, FLAG1 );
    tncg = tncg + ncg;
    
    %% Save iteration and continue

    
    % Measure the GN update length
    if count==1 
        
        rdm(count) =  1;
        gdm(1) = norm(m_in - invmod);
        
    else
        
        gdm(2) = norm(m_in - invmod);
        rdm(count) = abs( gdm(2) - gdm(1) ) / norm(invmod);
    
        gdm(1) = gdm(2);
        
    end
    
    phi_d(count) = sum(( G*invmod - d ).^2);
    phi_m(count) = (invmod)'*(MOF)*(invmod);
    phi(count) = objfunc(invmod,MOF,beta(count));
    
    if count ~= 1
        
        dphim(count) = abs(phi_m(count) - phi_m(count-1)) / phi_m(count) *100;
        
    end
    % Get truncated cells
    tcells = spdiags(Pac);
    
    

    fprintf(' phi_d:\t %8.5e \n',phi_d(count))
    fprintf(' phi_m:\t %8.5e \n',phi_m(count))
    fprintf(' dphi_m:\t %8.5e \n',dphim(count))

    
    fprintf('Number of Inactive cells: %i\n',sum(tcells));
    fprintf('Number of CGS iterations: %i\n\n',ncg);
    
end