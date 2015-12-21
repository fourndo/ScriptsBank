function [MOF,aVRWs,aVRWx,aVRWy,aVRWz] = get_lp_MOF_3D_v2(model,model_ref,phim,Ws,Wx,Wy,Wz,Gx,Gy,Gz,t,alpha,LP,FLAG1,FLAG2,FLAG3,delta_s,delta_xyz)
% Function get_lp_WRW(Ws,Wx,Wy,Wz,p,q,r)
% Generates the elements of the model objective function for a specified
% p-norm(m) and q-norm(grad m) and a ratio between the two (r)and for a
% model (m). The lp-norm is approximated through a lag diffusivity approach
% where:
% 
% V3: ALL SCALED ON PHI_l2 (2014-08-02) 
%
% gradphi = Ws' * Rs * Ws +  Wx' * Rx * Wx + ...
% R = diag(1./ ( f(m).^(2-p) + epsilon )
% f(m) is abs(grad(m)) or abs(m)
%
% INPUT:
% m  :   Model at the kth iteration
% Ws :   Weighting matrix for the smallest component
% Wx :   Weighting matrix for the gradient x-component
% Wy :   Weighting matrix for the gradient y-component
% Wz :   Weighting matrix for the gradient z-component
% delta = 1e-5;
mcell = size(Ws,1);

% Wr = spdiags(wr,0,mcell,mcell);

%% MODIFICATION FOR TMVI
% If LP is size(3*5), then the process is repeated for the three model
% spaces

nspace = size(LP,2)/5;



for jj = 1 : nspace
    
    % Temporary test - scale the phi and theta component for MVI
    if jj ==2
        
        scs = 1e-1;
        scx = 1e-1;
        
    elseif jj == 3
        
        scs = 1e-1;
        scx = 1e-1;
        
    else
        
        scs = 1;
        scx = 1;
        
    end
    
    d = zeros(1,nspace);
    d(jj) = 1;
    
    % Create sub-space matrix and grab the right model
    D = kron(d,speye(mcell));
    
    m = D*model;
    mref = D*model_ref;
    
    lp = LP(:,(1:5)+(jj-1)*5);
    
switch FLAG1
    case 'SMOOTH_MOD'
        % Compute gradients then compute magnitude         
        dmx = Gx * m;
        dmy = Gy * m;
        dmz = Gz * m;
        
    case 'SMOOTH_MOD_DIF'
        
        
        dmx = Gx * (m - mref);
        dmy = Gy * (m - mref);
        dmz = Gz * (m - mref);
        
    otherwise
        fprintf('FLAG must be SMOOTH_MOD | SMOOTH_MOD_DIF\n')
        fprintf('Please revise input\n')
        return
end


switch FLAG2
    case 'dmdx'
        % Compute gradients then compute magnitude 
        GRADx =  dmx;
        GRADy =  dmy;
        GRADz =  dmz;

    case 'GRADm'

        GRADx = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
        GRADy = sqrt(dmx.^2 + dmy.^2 + dmz.^2);
        GRADz = sqrt(dmx.^2 + dmy.^2 + dmz.^2);

    otherwise
        fprintf('FLAG must be GRADm | dmdx\n')
        fprintf('Please revise input\n')
        return
end

%% ADD SPECIAL SWITCH FOR MVI
% Impose l0 on the amplitude



%%

tRs = spalloc(mcell,mcell,mcell);
tRx = spalloc(mcell,mcell,mcell);
tRy = spalloc(mcell,mcell,mcell);
tRz = spalloc(mcell,mcell,mcell);




% if FLAG3~=0
    
    rvec = zeros(mcell,1);

    for ii = 1 : size(t,2)
        
        tt =t(:,ii);        
        
        tt(tt>0) = sqrt(tt(tt>0));
        
        T = spdiags((tt),0,mcell,mcell);
%         S = spdiags(s(:,ii),0,mcell,mcell);
        
%         if nspace > 1
% 
% 
%             mx = model(1:mcell);
%             my = model(1+mcell:2*mcell);
%             mz = model(1+2*mcell:3*mcell);
%             
%             mm = sqrt(mx.^2+my.^+mz.^2);
%             
%             rs = 1./ ( abs((mm - mref)) .^2 + delta_s.^2 ).^( (1 - lp(ii,1)/2) /2 );
%             
%         else    

            rs = 1./ ( abs((m - mref)) .^2 + delta_s(ii).^2 ).^( (1 - lp(ii,1)/2) /2 );

%         end

        
        rx = 1./ ( abs(GRADx) .^2 + delta_xyz.^2 ).^( (1 - lp(ii,2)/2) /2 );
        ry = 1./ ( abs(GRADy) .^2 + delta_xyz.^2 ).^( (1 - lp(ii,3)/2) /2 );
        rz = 1./ ( abs(GRADz) .^2 + delta_xyz.^2 ).^( (1 - lp(ii,4)/2) /2 );

        rvec = rvec + (T * ones(mcell,1)) * lp(ii,5);
        
        Rs = spdiags( rs ,0,length(rs),length(rs));

        Rx = spdiags( rx ,0,size(Wx,1),size(Wx,1));

        Ry = spdiags( ry ,0,size(Wy,1),size(Wy,1)); 
        Rz = spdiags( rz ,0,size(Wz,1),size(Wz,1));                  

%         subRs =  T * Rs ;
%         subRx =  T * Rx * Gx  ;
%         subRy =  T * Ry * Gy  ;
%         subRz =  T * Rz * Gz  ;
       
        
%         eta_s = 1 / max(abs((((subRs))'* ((subRs) * ((m - mref)) ) ) ) );
        eta_s = 1;%delta_s(ii).^(1 - lp(ii,1)/2);
        
%          switch FLAG1
%         case 'SMOOTH_MOD'
%             
%             % Compute gradients then compute magnitude 
% %             eta_x = 1 / max(abs( ( subRx'* ( (subRx *  m) ) ) ) );
% %             eta_y = 1 / max(abs( ( subRy'* ( (subRy *  m) ) ) ) );
% %             eta_z = 1 / max(abs( ( subRz'* ( (subRz *  m) ) ) ) );

        eta_x = 1;%delta_xyz(ii).^(1 - lp(ii,2)/2);
        eta_y = 1;%delta_xyz(ii).^(1 - lp(ii,3)/2);
        eta_z = 1;%delta_xyz(ii).^(1 - lp(ii,4)/2);
                
%         case 'SMOOTH_MOD_DIF'
% 
%             % Compute gradients then compute magnitude 
% %             eta_x = 1 / max(abs(  (subRx'* ( (subRx * ( m - mref ) ) ) ) ));
% %             eta_y = 1 / max(abs(  (subRy'* ( (subRy * ( m - mref ) ) ) ) ));
% %             eta_z = 1 / max(abs(  (subRz'* ( (subRz * ( m - mref ) ) ) ) ));
% 
%          end
        
        tRs = tRs + sqrt( eta_s ) * T * Rs ;
        tRx = tRx + sqrt( eta_x ) * T * Rx ;
        tRy = tRy + sqrt( eta_y ) * T * Ry ;
        tRz = tRz + sqrt( eta_z ) * T * Rz ;
        
        
             
    end
    
    avrws =  spdiags(sqrt( rvec * alpha(1) * scs ),0,mcell,mcell) * Ws * tRs ;
    avrwx =  spdiags(sqrt( (2 - rvec) * alpha(2) * scx ),0,mcell,mcell) * Wx * tRx *  Gx  ;
    avrwy =  spdiags(sqrt( (2 - rvec) * alpha(3) * scx ),0,mcell,mcell) * Wy * tRy *  Gy  ;
    avrwz =  spdiags(sqrt( (2 - rvec) * alpha(4) * scx ),0,mcell,mcell) * Wz * tRz *  Gz  ;                    
                        

    
% else
%        
%     % Form the matrices by including the weights and volume
%     avrws =  sqrt( alpha(1) ) * Ws ;
%     avrwx =  sqrt( alpha(2) ) * Wx * Gx;
%     avrwy =  sqrt( alpha(3) ) * Wy * Gy;
%     avrwz =  sqrt( alpha(4) ) * Wz * Gz;
%     
% end

% If more than one zone, average the scaling between lp,lq on transition
% if size(LP,1) ~= 1;
        
%     aVRWs = spdiags(sqrt( gamma * rvec ),0,mcell,mcell) * aVRWs;
%     aVRWx = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWx;
%     aVRWy = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWy;
%     aVRWz = spdiags(sqrt( gamma * ( 2.0 - rvec ) ) ,0,mcell,mcell) * aVRWz;
    
% Otherwise just a scaler
% else
%     
%     aVRWs = sqrt( gamma * LP(1,5) ) * aVRWs;
%     aVRWx = sqrt( sc_phigx * ( 2.0 - LP(1,5) ) ) * aVRWx;
%     aVRWy = sqrt( sc_phigy * ( 2.0 - LP(1,5) ) ) * aVRWy;
%     aVRWz = sqrt( sc_phigz * ( 2.0 - LP(1,5) ) ) * aVRWz;
%     
% end



    % Form the final model objective function
    if jj == 1
        
        aVRWs = kron(d,avrws);
        aVRWx = kron(d,avrwx);
        aVRWy = kron(d,avrwy);
        aVRWz = kron(d,avrwz);
        
    else
        
        aVRWs = [aVRWs;kron(d,avrws)];
        aVRWx = [aVRWx;kron(d,avrwx)];
        aVRWy = [aVRWy;kron(d,avrwy)];
        aVRWz = [aVRWz;kron(d,avrwz)];
        
    end

end

if FLAG3~=0
gamma = phim /(...
        ((aVRWs * (model-model_ref) )' * (aVRWs * (model-model_ref) ))+...
        (( (aVRWx) * (model) )' * ( (aVRWx) * (model) ))+...
        (( (aVRWy) * (model) )' * ( (aVRWy) * (model) ))+...
        (( (aVRWz) * (model) )' * ( (aVRWz) * (model) )));
% gamma=1;
aVRWs = sqrt(gamma) * aVRWs; 
aVRWx = sqrt(gamma) * aVRWx;
aVRWy = sqrt(gamma) * aVRWy;
aVRWz = sqrt(gamma) * aVRWz;
end

MOF = aVRWs'*aVRWs + aVRWx'*aVRWx + aVRWy'*aVRWy + aVRWz'*aVRWz;

