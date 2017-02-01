function [x,r,count]=CG_Lin(x,G, W,MOF,RHS, Proj, PreC, Patv)

r= Patv*RHS - (Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*x))) + Patv*MOF*Patv*x + Patv*x*1e-3 );    
p= PreC * r;

%Save sens weight to plot
% save(['C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1\SensW.dat'],'-ascii','p');




sold = r' * (PreC * r);
s0 = sold;
count=0;

% normr = norm(r);
while count < 20
    
    count=count+1;
       
    q= Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*p))) + Patv*MOF*Patv*p + Patv*p*1e-3;
    
    alpha=sold./(p'*q);
   
    x=x+alpha.*p;
     
    r=r-alpha*q; 
    
    h = PreC * r;
    snew=r'* h;
    p = h + (snew/sold.*p);
    
    sold=snew;

    if sqrt(snew / s0) < 1e-3
        
        fprintf('Number of CG %i:', count);
        return
        
    end
    
%     normr = norm(r);  

%     count=count+1;
    
%     figure(100)
%     semilogy (count,(snew) / (norm(s0)),'*');
%     hold on
%     
       
end

fprintf('Number of CG %i:', count);
% hold off