function [flag,expan,trans,zone,static,numcell,newdist]=fixtrans(flag,expan,trans,zone,static,numcell,newdist,counter)          

    while zone.trans(counter+1)>1.4 
        numcell=numcell+1;
        zone.trans(counter+1)=get_factor(numcell,zone.dz(counter+1),newdist);

        %Computes new transition factor
        trans.vec=zone.dz(counter+1)*zone.trans(counter+1).^(0:numcell-1);
    end

    while (zone.trans(counter+1)<1.0 || expan.vec(flag.expan)/trans.vec(end)>1.4) && numcell>=2
        if numcell==2;
            newdist=newdist+zone.dz(counter+1);      
            zone.trans(counter+1)=get_factor(numcell,zone.dz(counter+1),newdist);
            trans.vec=zone.dz(counter+1)*zone.trans(counter+1).^(1:numcell-1); 
            numcell=1;
        else
        numcell=numcell-1;

        zone.trans(counter+1)=get_factor(numcell,zone.dz(counter+1),newdist);
        trans.vec=zone.dz(counter+1)*zone.trans(counter+1).^(0:numcell-1);
        end
    end
