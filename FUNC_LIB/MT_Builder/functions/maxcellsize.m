function [zonemesh]=maxcellsize(trans,expan,zone,flag,node,static,jj)
% Computes the zone mesh if maximum cell size has been reached

% Test that max cell size is greater or equal to the cell size of the
% vectors




    %Find the number of cells in trans.vec and expan.vec greater
    %than the maxcell size
    countdown=0;

    while trans.vec(flag.trans-countdown)>zone.sizemax(jj) && flag.trans-countdown>1;
        countdown=countdown+1;
    end  

    countup=0;
    while expan.vec(flag.expan-countup)>zone.sizemax(jj) && flag.expan-countup>1;
        countup=countup+1;
    end 
    
% Compute the length taken by the vectors after zone.sizemax cells are
% removed
vectorlength=0;
for ii=1:(flag.trans-countdown)
vectorlength=vectorlength+trans.vec(ii);
end

for ii=1:(flag.expan-countup)
vectorlength=vectorlength+expan.vec(ii);
end

%Number of zone.sizemax cells
 numcellmax=floor((zone.depth(jj)- vectorlength) / zone.sizemax(jj));

 if numcellmax==0
     numcellmax=1;
 end
 
% Start building the new mesh using the expansion vector and the
% zone.sizemax. The transition vector will take priority if it the cell 
% size is equal to sizecell max
    zonemesh=[expan.vec(1:(flag.expan-countup)) ones(1,numcellmax)*zone.sizemax(jj)];
    
%Compute a new transition vector to complete the zone
numcell=flag.trans-countdown;
if numcell<=1
    numcellmax=numcellmax-1;
    zonemesh=[expan.vec(1:(flag.expan-countup)) ones(1,numcellmax)*zone.sizemax(jj)];
    numcell=numcell+1;
end
newdist=zone.depth(jj)-sum(zonemesh);

zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);

% Check for good transition
while zone.trans(jj+1)>1.4 
    numcell=numcell+1;
    zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);

    %Computes new transition factor
    trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
end

while (zone.trans(jj+1)<1.0 || zone.sizemax(jj)/trans.vec(end)>1.4) && numcell>=2
    if numcell==2;
        newdist=newdist+zone.dz(jj+1);      
        zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
        trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(1:numcell-1); 
        numcell=1;
    else
    numcell=numcell-1;

    zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
    trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
    end
end

% If the transition is still too step, then the expan flag is moved down of one.   
while (zone.trans(jj+1)>1.4 || zone.trans(jj+1)<1.0 || zone.sizemax(jj)/trans.vec(end)>1.4) && numcellmax>=1
    numcellmax=numcellmax-1;
    zonemesh=[expan.vec(1:(flag.expan-countup)) ones(1,numcellmax)*zone.sizemax(jj)];
    
    newdist=zone.depth(jj)-sum(zonemesh);
    numcell=flag.trans;

    zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
    trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);

    while zone.trans(jj+1)>1.4 
        numcell=numcell+1;
        zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);

        %Computes new transition factor
        trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
    end

    while (zone.trans(jj+1)<1.1 || zone.sizemax(jj)/trans.vec(end)>1.4) && numcell>=2
        if numcell==2;
            newdist=newdist+zone.dz(jj+1);      
            zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
            trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(1:numcell-1); 
            numcell=1;
        else
        numcell=numcell-1;

        zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
        trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
        end
    end
end


%     %Zone mesh
    zonemesh=[expan.vec(1:(flag.expan-countup)) ones(1,numcellmax)*zone.sizemax(jj) trans.vec(end:-1:1)];