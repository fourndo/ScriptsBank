function [zone,static,node,code]=zonecheck(zone,static,node)

code=0;


%% Test #1  - Remove small zones

count=2;
while count<length(zone.dz)-1;
    
    % Find smallest cell size between zone above and below. The small zone
    % is merged to the zone with smallest cell size.
    tempo = find(([zone.dz(count-1) zone.dz(count+1)])==min([zone.dz(count-1) zone.dz(count+1)]));
    
    % Keep the top zone in case if both zones above and below have the same
    % cell size. Flag will have value of -1 or 1 (above or below).
    
    flag = tempo(1); %If above and below have same cell size, select zone above.
    
    if flag == 1
        flag = -1;
    else
        flag = 1;
    end
    
    if zone.depth(count) < 5 * zone.dz(count) %&& (zone.rank(count)<zone.rank(count+flag))
       % Skrink the parameters to get rid of small zone and move nodes and
       % zone.depth to adapt the model.
       zone.depth(count+flag) = zone.depth(count+flag) + zone.depth(count);
       zone.depth = zone.depth([1:(count-1) (count+1):end]);
       zone.dz = zone.dz([1:(count-1) (count+1):end]);
       zone.expan = zone.expan([1:(count-1) (count+1):end]);
       zone.trans = zone.trans([1:(count-1) (count+1):end]);
       zone.sizemax = zone.sizemax([1:(count-1) (count+1):end]);
       
       %Keep the right node!!
        if flag == -1
            node = node([1:(count-1) (count+1):end]);
        else
            node = node([1:(count) (count+2):end]);
        end
       
       sprintf('Zone %i has been merged to zone %i',count,count+flag)
       
    end
    count=count+1;
end

%% Test #2  - Push boundary of zone down if expansion = 1.0 to get interger number of
% cells
for ii=1:length(zone.dz)-1;
    
    
    if zone.expan(ii) == 1.0
       
       depth_tempo = ceil(zone.depth(ii)/zone.dz(ii))*zone.dz(ii);
       zone.depth(ii+1) = zone.depth(ii+1) - (depth_tempo-zone.depth(ii));
       zone.depth(ii) = depth_tempo;
       node(ii+1) = node(ii+1) - (depth_tempo-zone.depth(ii));
       sprintf('Zone depth %i has been changed to %4.1f to finish exactly',ii,zone.depth(ii))
       
    end
    
end

%% Test #3  - Sanity check on cell sizes for easy transition
for ii=2:length(node)-2
%     if ii==1 && zone.depth(ii)/zone.dz(ii+1)<3
% %         sprintf('Zone %i has only %4.1f cells, need at least 3',ii,zone.depth(ii)/zone.dz(ii+1))
% %         sprintf('Cell size of zone %i has been changed to %4.1f m', ii, (zone.depth(ii)/3))
% %         zone.dz(ii)=(zone.depth(ii)/3);
%         
%     elseif ii==length(node)-1 && zone.depth(ii)/zone.dz(ii-1)<3
% %         sprintf('Zone %i has less than one cell from above',ii,zone.depth(ii)/zone.dz(ii))
% %         sprintf('Cell size of zone %i has been changed to %4.1f m', ii, (zone.depth(ii)/3))
% %         zone.dz(ii-1)=(zone.depth(ii)/3);
        
    if zone.depth(ii)/zone.dz(ii)<5;
        sprintf('Zone %i has only %4.1f cells, need at least 5',ii,zone.depth(ii)/zone.dz(ii))
        sprintf('Cell size of zone %i has been changed to %4.1f m', ii, (zone.depth(ii)/5))
        zone.dz(ii)=(zone.depth(ii)/5);
    end
end

%% Test #4  - Test on the maxcellsize must be at least as big as the cell expanding
% from.
for ii=1:length(zone.sizemax)
    if ii==1 && zone.sizemax(ii)<zone.dz(ii+1)
        zone.sizemax(ii) = zone.dz(ii+1) * zone.expan(ii);
        sprintf('Zone %i has maximum cell size too small',ii)
        sprintf('Maximum cell size of zone %i has been changed to %4.1f m', ii, zone.sizemax(ii)) 
    
    elseif ii==length(zone.sizemax) && zone.sizemax(ii)<zone.dz(ii-1)
        zone.sizemax(ii) = zone.dz(ii-1) * zone.expan(ii);
        sprintf('Zone %i has maximum cell size too small',ii)
        sprintf('Maximum cell size of zone %i has been changed to %4.1f m', ii, zone.sizemax(ii))
        
    elseif ii~=length(zone.sizemax) && ii~=1 && zone.sizemax(ii) < zone.dz(ii)
        zone.sizemax(ii) = zone.dz(ii) * zone.expan(ii) ;
        sprintf('Zone %i has maximum cell size too small',ii)
        sprintf('Maximum cell size of zone %i has been changed to %4.1f m', ii, zone.sizemax(ii))      
    end
end

% %% Test #5  -  Check if first and last zone have enough space to transition
% for ii = [1 length(node)-1]    
%     if ii==1 && zone.depth(ii) < sum( zone.dz(ii+1) * zone.trans(ii) .^(0:2) )
%        sprintf('Zone %i is to shallow of transitions, please verify',ii)
%        code=1;
%        break
%     end
% 
%     if ii==length(node)-1 && zone.depth(ii) < sum( zone.dz(ii-1) * zone.trans(ii) .^(0:2) )
%        sprintf('Zone %i is to shallow of transitions, please verify',ii)
%        code=1;
%        break
%     end
% end

% %% Test #6  -  Check if we can transition from below leaving at least 3 cells
% % for all zones except the top and bottom
% for ii=2:length(node)-2
% 
%     
%     if zone.dz(ii)>=zone.dz(ii+1)
%         
%         newdist=( zone.depth(ii) - 3 * zone.dz(ii));
%         trans.length=0;
%         counter=0;
%         
%             if newdist<0
%                sprintf('Zone %i is to shallow of transitions, please verify',ii)
%                code=1;
%                break
%             end
%             
%             while trans.length < newdist
%             counter=counter+1;
%             trans.vec(counter)=zone.dz(ii+1)*zone.trans(ii)^(counter-1);
%             trans.length=trans.length+trans.vec(counter);
%             end
% 
%         step=zone.dz(ii)/trans.vec(counter-1);
%         
%         if (step)>1.4;
%         sprintf('Transition between zone %i and %i is too steep',ii,ii+1)      
%         %Compute new cell size for zone ii
%         zone.dz(ii)=trans.vec(counter-1);    
%         sprintf('Cell size of zone %i has been changed to %4.1f m', ii, zone.dz(ii)) 
%         end
%         
%     
%     
% 
%     
%     elseif zone.dz(ii)<=zone.dz(ii+1) 
%         
%         newdist=(zone.depth(ii+1)-3*zone.dz(ii+1)); %Keep at least three cells
%         trans.length=0;
%         counter=0;
%             if newdist<0
%                sprintf('Zone %i is to shallow of transitions, please verify',ii+1)
%                code=1;
%                break
%             end
%             
%             while trans.length < newdist
%             counter=counter+1;
%             trans.vec(counter)=zone.dz(ii)*zone.trans(ii+1)^(counter-1);
%             trans.length=trans.length+trans.vec(counter);
%             end
% 
%         step=zone.dz(ii+1)/trans.vec(counter);
%         
%         
%         if(step)>1.4;
%         sprintf('Transition between zone %i and %i is too steep',ii,ii+1)      
%         %Compute new cell size for zone ii
%         zone.dz(ii+1)=trans.vec(counter-1);    
%         sprintf('Cell size of zone %i has been changed to %4.1f m', ii, zone.dz(ii))
%         end
%         
%     end
%     
% end

