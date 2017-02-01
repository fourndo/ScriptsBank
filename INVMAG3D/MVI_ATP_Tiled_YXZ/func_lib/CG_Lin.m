function [x,r,count]=CG_Lin(x,G, W,MOF,RHS, Proj, PreC, Patv)

r= Patv*RHS - (Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*x))) + Patv*MOF*Patv*x);    
p= PreC * r;

%Save sens weight to plot
% save(['C:\Users\dominiquef.MIRAGEOSCIENCE\ownCloud\Research\Modelling\Synthetic\Nut_Cracker\Tiled_CMI\Tile1\SensW.dat'],'-ascii','p');


s0 = p;

sold = r' * s0;
s0=sold;
count=0;

% normr = norm(r);
while count < 200
    
    count=count+1;
       
    q= Patv*Proj'*Gtvec(G,W,Gvec(G,W,Proj*(Patv*p))) + Patv*MOF*Patv*p;
    
    alpha=sold./(p'*q);
   
    x=x+alpha.*p;
     
    r=r-alpha*q; 
    
    h = PreC * r;
    snew=r'* h;
    p = h + (snew/sold.*p);
    
    sold=snew;

    if (snew) / ((s0)) < 1e-3
        fprintf('CG steps %i: \n',count);
        return
        
    end
    
%     normr = norm(r);  

%     count=count+1;
    
%     figure(100)
%     semilogy (count,(snew) / (norm(s0)),'*');
%     hold on
%     
       
end
fprintf('CG steps %i: \n',count);
% hold off