function [flag,expan,trans,zone,static,numcell,newdist]=fixexpan(flag,expan,trans,zone,static,numcell,newdist,counter)          



while zone.expan(counter)>1.4 || expan.vec(end) > zone.sizemax(counter)
    numcell=numcell+1;
    zone.expan(counter)=get_factor(numcell,zone.dz(counter),newdist);
    expan.vec=zone.dz(counter)*zone.expan(counter).^(0:numcell-1);
end

while (zone.expan(counter)<1.0 || zone.sizemax(counter)/expan.vec(end)>1.4) && numcell>=2

    %Get rid of static cells if not enough space
    while flag.expan<=2 && static(counter)~=0
        numcell= flag.expan + 1;
        static(jj)= static(jj) - 1;
        newdist= zone.depth(jj) - flag.trans * zone.sizemax(jj) - static(jj) * zone.dz(jj);
    end

    if numcell==2;
        newdist=newdist+zone.dz(counter);      
        zone.expan(counter)=get_factor(numcell,zone.dz(counter),newdist);
        expan.vec=zone.dz(counter)*zone.expan(counter).^(1:numcell-1); 
        numcell=1;
    else

    numcell=numcell-1;

    zone.expan(counter)=get_factor(numcell,zone.dz(counter),newdist);
    expan.vec=zone.dz(counter)*zone.expan(counter).^(0:numcell-1);
    end
end