function [m,tncg] = GNsolver( G, m, d, phi, MOF, beta , aVRWs, aVRWx, aVRWy, aVRWz )
% Solve H dm = J and update model

mcell = length(m);

tncg = 0;        % Count total CG steps
solves = 1;       % Count number of GN solves
rddm = 1;       % Mesure the relative change in rel|dm|
ddm = [1 1];    % Mesure of change in model update |dm|
 
objfunc = @(x) sum( ( G * x - d ).^2 ) + ( x' * beta * MOF * x );
cg_count = 0;
while solves < 10 && rddm > 1e-4

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
                [dm,~,ncg] = CGLSQ( dm, A , b );
                
                tncg = tncg + ncg;
                %% Step length, line search                
%                 ncg = ncg+iter; % Record the number of CG iterations

                gamma = 2;

                % Reduce step length in order to reduce phid
                phi_temp = 0;   
                while (phi_temp > phi || gamma == 2) && gamma>1e-4

                    gamma = 0.5 * gamma;

                    gdm = gamma * dm;

                    ddm(2) = norm(dm);

                    m_temp = m + gdm;

                    phi_temp = objfunc(m_temp);

                end
              
                % Mesure relative change in dm
                if solves == 1

                    rddm = 1;
                    ddm(1) = ddm(2);

                else

                    rddm = ddm(2)/ddm(1);

                end

                
                % Update model
                m = m + dm;
                
%                 fprintf('GN iter %i |g| rel:\t\t %8.5e\n',solves,rddm);
%                 fprintf('CG iter %i \n',ncg);
                solves = solves + 1;

end                