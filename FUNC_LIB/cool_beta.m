function [arg_out,beta_out] = cool_beta( beta_in , phi_d , dphim, target , arg_in , tol, dphim_min)
% Function cool beta for different stages of inversion
% INPUT: 
% beta_in : Current beta
% targer: target phi_d
% phi_d: current misfit
% tol: tolerance in percent
% 
% Ouput
% arg: For different stages
% beta_out

% If model update is smaller than tolerance, then finish the inversion
if (dphim < dphim_min && arg_in == 2) && (phi_d > target * (1-tol) && phi_d < target * (1+tol) && arg_in == 2)

    arg_out = 3;
    beta_out = beta_in;
    fprintf('# FINAL ITERATION #\n');

% Otherwise change beta
else

    % Check to see if overshooted the target misfit,
    % If yes, then switch to mode 1 and redo the iteration as an lp-lq
    % Mode 1: Lp-lq inversion with smaller beta steps
    if phi_d < target * (1-tol) %&& arg_in~=1                    

        

        if arg_in == 0
            fprintf('# NEXT Lp-Lq - STATIONNARY BETA #\n');
%             beta_out = beta_in / (1 - tol) ;
            beta_out = beta_in;            
            arg_out = 1;

        else

fprintf('# NEXT ITER - INCREASING BETA #\n');
             beta_out = beta_in  * target / phi_d ;
%             beta_out = beta_in / (1 - tol) ;
            arg_out = arg_in;

        end


    % Misfit is still too large -> reduce beta
    elseif phi_d > target * (1 + tol) %&& arg_in ~= 1


        fprintf('# NEXT ITER - REDUCING BETA #\n');
        if arg_in == 0

            beta_out = 0.75 * beta_in; 
            arg_out = arg_in;

        else

            beta_out = beta_in  * target / phi_d ;
            arg_out = arg_in;

        end               


    else % Keep same beta and iterate again

        fprintf('# NEXT ITER - STATIONNARY STEP #\n');

        if arg_in == 0

            arg_out = 1;

        else

            arg_out = arg_in;

        end

        beta_out = beta_in;


    end
    
end
