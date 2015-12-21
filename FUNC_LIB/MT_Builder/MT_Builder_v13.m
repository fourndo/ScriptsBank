% Mira Geoscience
% Project 3802_FugroEM - BETA VERSION v12.0
% Adaptative mesh grid with transition ans expansion factors
% Written by: D Fournier
% Date: January 4th, 2012


clear all
close all

addpath functions
%% Step 1 - User Inputs

% For this specific example, the inputs are:
zmax=0;
zmin=-2000;
node=[zmax zmin];


for ii=1:length(node)-1
zone.trans(ii)=1.4;                 %Transition factors out of the zones
% zone.expan(ii)=1.4;                 %Expansion factors within each zones
zone.depth(ii)=node(ii)-node(ii+1); %Zone length vertical
% sprintf('Maximum cell size of %3.2f in zone %i',zone.depth(ii)/5,ii)
zone.sizemax(ii)=100;
static(ii)=0;
end

zone.expan(1)=1.2;
zone.expan(2)=1.1;
zone.expan(3)=1.2;
% zone.expan(4)=1.0;
% zone.expan(5)=1.05;
% zone.expan(6)=1.0;
% zone.expan(7)=1.2;


%User specify maximum cell size per zone
zone.dz(1)=100;
zone.dz(2)=12.5;
zone.dz(3)=69.499;
% zone.dz(4)=75;
% zone.dz(5)=15;
% zone.dz(6)=15;
% zone.dz(7)=25;
% zone.dz(8)=25;

%User specify maxcellsize
zone.sizemax(1)=50;
zone.sizemax(2)=100;
zone.sizemax(3)=83.3988;
% zone.sizemax(4)=25;
% zone.sizemax(5)=15;
% zone.sizemax(6)=400;


 % Number of cells that are not expanded
code=0;

%User set zone priority
zone.rank(1)=0;
zone.rank(2)=0;
zone.rank(3)=0;
zone.rank(4)=0;
zone.rank(5)=0;
zone.rank(6)=0;

% Check that transitions are possible between zones and that the chosen
% cell size is feasible.
if length(node)>3
[zone,static,node,code]=zonecheck(zone,static,node);
end

if code==1
    break
end

% Initialising
mesh=[];
trackzone=[];
jj=1;

% Iterate over all zones
while jj <= length(node)- 1
    
    
%% Compute first and last zone
% First zone is computed differently than the rest as we use the 
% zone.sizemax as the transition vector and the cell size of the zone
% below for the expansion. If there is only one zone, the zonemesh is
% flipped for the zone.size max to be at the bottom.

if jj == 1 || ( jj== length(node) - 1 && length(node)>2 )
    
    if jj==1
       zone.dz(jj)= zone.dz(jj+1);
    else
       zone.dz(jj)= mesh(end); 
    end
    
    expan.vec= zone.dz(jj);  %Initialize expansion vector cells
    expan.length= expan.vec;   %Initialize expansion vector length
    trans.vec= zone.sizemax(jj);
    trans.length= trans.vec;
    
    expan.marker= expan.vec / 2;   %Markers for expansion vector cells
    trans.marker= zone.depth(jj) - trans.vec / 2;   %Markers for transition vector cells
    
    count=1;
        
 	%Build expan and trans vectors                             
        while expan.length(count)<zone.depth(jj) || trans.length(count)<zone.depth(jj)
        count=count+1;
        
        if count<=static(jj)
        expan.vec(count)=zone.dz(jj);
        else
        expan.vec(count)= zone.dz(jj) * zone.expan(jj) ^ (count - 1 - static(jj));
        end
        trans.vec(count)= zone.sizemax(jj);
        
        expan.length(count)= sum(expan.vec);
        trans.length(count)= sum(trans.vec);
        
        expan.marker(count)= expan.length(count-1) + expan.vec(count) / 2;
        trans.marker(count)= zone.depth(jj) - trans.length(count-1) - trans.vec(count) / 2;       
        end
    
    % Find best overlap
    [flag.expan,flag.trans]= get_overlap(expan.vec,trans.vec,expan.marker,trans.marker);
    
    % Deal with first zone that has an expansion factor of 1.0
     if zone.expan(jj) == 1.0
        numcell= flag.expan;
        % Expansion vector takes the left-over from the zone
        remainder= zone.depth(jj) - expan.length(numcell);
            
            while remainder > expan.vec(numcell)
                numcell= numcell + 1;
                expan.vec(numcell)= expan.vec(numcell-1) * zone.expan(jj);
                expan.length(numcell)= sum(expan.vec(1:numcell));
                remainder= zone.depth(jj) - expan.length(numcell);
            end
            
            if remainder<0
                numcell= numcell - 1;
              remainder= zone.depth(jj) - expan.length(numcell); 
            end
            
            %Absorb the remainder in the cell then distribute the remainder
            %linearly through the cells to reduce the jump
%             expan.vec=[expan.vec(1:numcell-1) expan.vec(numcell)+remainder];
            
            buffercells= 1;
            buffer= expan.vec(numcell) + remainder;
            while buffer(1) / expan.vec(end-buffercells)>1.4 && (numcell-buffercells)>1
                buffercells= buffercells + 1;
                
                buffer= expan.vec(numcell-buffercells+1:numcell) +...
                    [1:buffercells] .*remainder / ( sum (1:buffercells) );
                
            end
            expan.vec= [expan.vec(1:(numcell-buffercells)) buffer];
            
            % If there is only one zone, the zonemesh is flipped
            if length(node)==2 || jj==length(node)-1
            zonemesh=expan.vec;
            else
            zonemesh=expan.vec(end:-1:1); 
            end
            
    % If the expansion has an expansion factor but doesn't reach the 
    % zone.maxsize, the expansion fills the zone entirely.        
    elseif flag.trans==1; 
        
        % Case where the zone has too few space, then the remainder becomes
        % a cell. Might want to find something better, but easy fix for
        % now.
        if flag.expan==1;
            
           zonemesh = zone.depth(jj);
           sprintf('Zone %i is to shallow of transitions, please verify',jj)
           
        else

        numcell=flag.expan;
        newdist=zone.depth(jj)-static(jj)*zone.dz(jj);
        
        

        zone.expan(jj)= get_factor(numcell,zone.dz(jj),newdist);
        expan.vec= zone.dz(jj) * zone.expan(jj) .^(0:numcell-1);
        
        zone.sizemax(jj)=expan.vec(end);
        
        [flag,expan,trans,zone,static,numcell,newdist]=fixexpan...
            (flag,expan,trans,zone,static,numcell,newdist,jj);
        
        while zone.expan(jj)>1.4 && static(jj)~=0
            %Get rid of static cells if not enough space
            while numcell<=2 && static(jj)~=0
                numcell= numcell + 1;
                static(jj)= static(jj) - 1;
                newdist=zone.depth(jj)-static(jj)*zone.dz(jj);
            end

            zone.expan(jj)=get_factor(numcell,zone.dz(jj),newdist);
            expan.vec=zone.dz(jj)*zone.expan(jj).^(0:numcell-1);


            [flag,expan,trans,zone,static,numcell,newdist]=fixexpan...
            (flag,expan,trans,zone,static,numcell,newdist,jj);
        end
        % If there is only one zone, the zonemesh is flipped
            if length(node)==2 || jj==length(node)-1
            zonemesh=[ones(1,static(jj))*zone.dz(jj) expan.vec];
            else
            zonemesh=[expan.vec(end:-1:1) ones(1,static(jj))*zone.dz(jj)]; 
            end
    
        end
    % Otherwise the transition keeps the size zone.maxsize cell and
    % the expansion adapts
    else
        
    trans.vec=ones(1,flag.trans) * zone.sizemax(jj);
    newdist= zone.depth(jj) - flag.trans * zone.sizemax(jj) - static(jj) * zone.dz(jj);
    numcell= flag.expan - static(jj);
 
    %Get rid of static cells if not enough space
    while flag.expan<=2 && static(jj)~=0
        numcell= flag.expan + 1;
        static(jj)= static(jj) - 1;
        newdist= zone.depth(jj) - flag.trans * zone.sizemax(jj) - static(jj) * zone.dz(jj);
    end
    
    
    % Deal with case where best overlap is at the first cell of expansion
    % vector. Move all the flags and recompute
    while flag.expan<=1 && flag.trans>1
        flag.expan= flag.expan + 1;
        flag.trans= flag.trans - 1;
        numcell= flag.expan;
        trans.vec=ones(1,flag.trans) * zone.sizemax(jj);
        newdist= zone.depth(jj) - flag.trans * zone.sizemax(jj) - static(jj) * zone.dz(jj);
    end
    
        
    zone.expan(jj)=get_factor(numcell,zone.dz(jj),newdist);
    expan.vec=zone.dz(jj)*zone.expan(jj).^(0:numcell-1);
    

    [flag,expan,trans,zone,static,numcell,newdist]=fixexpan...
        (flag,expan,trans,zone,static,numcell,newdist,jj);

% If the transition is still too step, then we removed one sizemax cell
% until the transition gets smooth.
     while (zone.expan(jj)>1.4 || zone.expan(jj)<1.0 ||...
             zone.sizemax(jj)/expan.vec(end)>1.4 ||...
             expan.vec(end) > zone.sizemax(jj)) && flag.trans>1
         
        flag.trans=flag.trans-1;
        newdist= zone.depth(jj) - flag.trans*zone.sizemax(jj) -...
            static(jj) * zone.dz(jj);
        
        numcell=flag.expan - static(jj);

        zone.expan(jj)=get_factor(numcell,zone.dz(jj),newdist);
        expan.vec=zone.dz(jj)*zone.expan(jj).^(0:numcell-1);               

        [flag,expan,trans,zone,static,numcell,newdist]=fixexpan...
            (flag,expan,trans,zone,static,numcell,newdist,jj);

     end
            
            % If there is only one zone or it is the bottom zone, the zonemesh is flipped
            if length(node)==2 || jj==length(node)-1
            zonemesh=[ones(1,static(jj))*zone.dz(jj) expan.vec ones(1,flag.trans)*zone.sizemax(jj)];
            else
            zonemesh=[ones(1,flag.trans)*zone.sizemax(jj) expan.vec(end:-1:1) ones(1,static(jj))*zone.dz(jj)];
            end
    
    
    end


    mesh=[mesh zonemesh];
    trackzone=[trackzone ones(1,length(zonemesh))*jj];
    
     
    
else  
%% Step 2 - Iterate over all other zones
% Following program cycles through each zone and compute the expansion
% vector for the zone and the transition vector from the zone below upward.


    code=1;
    expan.vec= zone.dz(jj);  %Initialize expansion vector cells
    trans.vec= zone.dz(jj+1);%Initialize transition vector cells
    
    expan.length=expan.vec;   %Initialize expansion vector length
    trans.length=trans.vec;   %Initialize transition vector length
    
    expan.marker=expan.vec/2;   %Markers for expansion vector cells
    trans.marker=zone.depth(jj)-trans.vec/2;   %Markers for transition vector cells
    
%     if zone.sizemax(jj) <  trans.vec(1)*zone.trans(jj+1);
%     zone.sizemax(jj) = trans.vec(1)*zone.trans(jj+1);
%     elseif zone.sizemax(jj) < expan.vec(1)*zone.expan(jj);
%         zone.sizemax(jj) = expan.vec(1)*zone.expan(jj);
%     end
    
    %Compute expansion vector and transition from below through the whole zone.
    count=1;
    while expan.length(count)<zone.depth(jj) || trans.length(count)<zone.depth(jj)
        count=count+1;
        
        if count<=static(jj)
        expan.vec(count)=zone.dz(jj);
        else
        expan.vec(count)=zone.dz(jj)*zone.expan(jj)^(count-1-static(jj));
        end
        trans.vec(count)=zone.dz(jj+1)*zone.trans(jj+1)^(count-1);
        
        expan.length(count)=sum(expan.vec);
        trans.length(count)=sum(trans.vec);
        
        expan.marker(count)=expan.length(count-1)+expan.vec(count)/2;
        trans.marker(count)=zone.depth(jj)-trans.length(count-1)-trans.vec(count)/2;
        
    end

    
%% Case 1 - Zone below has a smaller cell size
% The expansion vector keeps its length and the transition vector adapts
    if zone.dz(jj+1)<zone.dz(jj)     
        
    [flag.expan,flag.trans]=get_overlap(expan.vec,trans.vec,expan.marker,trans.marker);
    
        %If the best overlapping cell is larger than the maxcellsize, than
        %the middle cells are replaced by maxcellsize and the vectors adapt
        %to it.
        if expan.vec(flag.expan)>zone.sizemax(jj) || trans.vec(flag.trans)>zone.sizemax(jj)

            [zonemesh,zone,node]=maxcellsize(trans,expan,zone,flag,node,static(jj),jj);
            mesh=[mesh zonemesh];
            trackzone=[trackzone ones(1,length(zonemesh))*jj];
            jj=jj+1;
            continue
        else

        %Compute new distance to fill by trans.vec and the nb of cells
        %to reach it
        newdist= zone.depth(jj) - expan.length(flag.expan);
        numcell= flag.trans;

        %If the transition is too long, make the large cells
        %constant and recompute the get_overlapping cell. A minimum of
        %three cells must be left to expansion to allow transition from
        %the upper zone.
        if flag.expan<=3 
            expan.vec(:)= zone.dz(jj);
            expan.marker= ( (0:length(expan.vec) - 1) + 0.5 ) * zone.dz(jj);

            %Recompute best overlap
            [flag.trans,flag.expan]=get_overlap(trans.vec,expan.vec,trans.marker,expan.marker);
%                 
            numcell= flag.trans;
            newdist= zone.depth(jj) - flag.expan*zone.dz(jj);
        end
        
                    
        % If the best overlap does not give enough space for the transition
        % vector, the overlap is moved up
        if flag.trans <=2 && flag.expan > 3          
            newdist = newdist + expan.vec(flag.expan);
            numcell = numcell + 1;
            flag.expan = flag.expan - 1;
            flag.trans = flag.trans + 1;
        end
        
% Back to General Case 1                 
            
        zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
        trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
        
        
        
        % Check if transition is between 1.0 and 1.4   
        [flag,expan,trans,zone,static,numcell,newdist]=fixtrans...
            (flag,expan,trans,zone,static,numcell,newdist,jj);
                
               
        % If the transition is still too step, then the expan flag
        % is moved down of one until smooth transition. 
        while (zone.trans(jj+1)>1.4 || zone.trans(jj+1)<1.0 ||...
                expan.vec(flag.expan)/trans.vec(end)>1.4) && flag.expan>2
            
            flag.expan=flag.expan-1;
            newdist=zone.depth(jj)-expan.length(flag.expan);
            numcell=flag.trans+1;

            zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
            trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);
            
            [flag,expan,trans,zone,static,numcell,newdist]=fixtrans...
                (flag,expan,trans,zone,static,numcell,newdist,jj);
        end
                
        end
        zonemesh=[expan.vec(1:flag.expan) trans.vec(end:-1:1)];


        mesh=[mesh zonemesh];
        trackzone=[trackzone ones(1,length(zonemesh))*jj];
%         sprintf('Zone %i starts at %6.3f',jj+1,zmax-sum(mesh))
        

    else

%% Case 2 - Zone below has a larger cell size
% The transition occurs within the current zone if and only if the
% expansion vector allows it. Otherwise, the expansion completes the zone,
% and the transition occurs in the zone below

    [flag.trans,flag.expan]=get_overlap(trans.vec,expan.vec,trans.marker,expan.marker);

    newdist=zone.depth(jj)-expan.length(flag.expan);
    numcell=flag.trans;
    

% If last expansion cell is smaller than the cells in the zone below.
% The program expands the mesh into the lower zone and moves the
% boundary down. This is generally the case for an expansion factor near 1.0    
        if flag.trans==1
            numcell=flag.expan;
            % Expansion vector takes the left-over from the zone
            remainder=zone.depth(jj)-expan.length(numcell);
            
            while remainder>expan.vec(numcell)
                numcell=numcell+1;
                expan.vec(numcell)=expan.vec(numcell-1)*zone.expan(jj);
                expan.length(numcell)=sum(expan.vec(1:numcell));
                remainder=zone.depth(jj)-expan.length(numcell);
            end
            
            if remainder<0
                numcell=numcell-1;
                remainder=zone.depth(jj)-expan.length(numcell); 
            end
            
            %Absorb the remainder in the cell then distribute the remainder
            %linearly through the cells to reduce the jump
%             expan.vec=[expan.vec(1:numcell-1) expan.vec(numcell)+remainder];
            
            buffercells=0;
            buffer=remainder;
            
            while (buffer(1)/expan.vec(numcell-buffercells)>1.4 ||...
                    buffer(1)/expan.vec(numcell-buffercells)<(1/1.4)) &&...
                    (numcell-buffercells)>1
                
                buffercells=buffercells+1;
                
                if buffercells==1;
                buffer=expan.vec(numcell)+remainder;
                else
                    
                buffer=expan.vec(numcell-buffercells+1:numcell)+[1:buffercells].*remainder/(sum(1:buffercells));
                end
                
            end
            expan.vec=[expan.vec(1:(numcell-buffercells)) buffer];
            
            
            % Generate a new transition vector from the end of the
            % expansion factor.    
        if expan.vec(end)<=zone.dz(jj+1) && zone.dz(jj+1)/expan.vec(end)>1.4
            transtart= expan.vec(end) * zone.trans(jj+1);
            trans.vec= transtart;
            trans.length = trans.vec;
            count=1;
            while trans.vec(count)<=zone.dz(jj+1) && trans.length(count)<(zone.depth(jj+1)-3*zone.dz(jj+1)) && zone.dz(jj+1)/trans.vec(count)>1.4
                count=count+1;
                trans.vec(count)=trans.vec(count-1)*zone.trans(jj+1);
                trans.length(count)=sum(trans.vec);
                
            end
            
            %Compute number of cells that can be taken away and recompute
            %the transition vector to match the distance exactly
            takecell=ceil(trans.length(count)/zone.dz(jj+1));
            newdist=takecell*zone.dz(jj+1);
            numcell=count;
            
            if numcell==1
              trans.vec=newdist;
            else
            zone.trans(jj+1)=get_factor(numcell,trans.vec(1),newdist);
            trans.vec=trans.vec(1)*zone.trans(jj+1).^(0:numcell-1);
            end
            
            % Check to make sure the transition has not exceeded the zone
            % below
            if numcell==1 && newdist > zone.depth(jj+1)
               
                trans.vec = zone.depth(jj+1);
                
            else
            %Check transition
            while (zone.trans(jj+1)>1.4 || trans.vec(end) > zone.dz(jj+1) || trans.vec(1)/expan.vec(end)>1.4) && zone.trans(jj+1)>1.0
                numcell=numcell+1;
                zone.trans(jj+1)= get_factor(numcell,transtart,newdist);

                %Computes new transition factor
                trans.vec= transtart*zone.trans(jj+1).^(0:numcell-1);
            end
            
            while (zone.trans(jj+1)<1.0 || zone.dz(jj+1)/trans.vec(end)>1.4) && numcell > 2
                numcell=numcell-1;
                zone.trans(jj+1)=get_factor(numcell,transtart,newdist);
                trans.vec= transtart * zone.trans(jj + 1).^(0:numcell-1);               
            end

            while (zone.trans(jj+1)<1.0 || zone.trans(jj+1)>1.4 || zone.dz(jj+1)/trans.vec(end)>1.4 || trans.vec(end) > zone.dz(jj+1)) && newdist <(zone.depth(jj+1)-3*zone.dz(jj+1))
                if newdist <(zone.depth(jj+1)-3*zone.dz(jj+1));
                    newdist=newdist+zone.dz(jj+1);      
                    zone.trans(jj+1)=get_factor(numcell,transtart,newdist);
                    trans.vec=transtart*zone.trans(jj+1).^(0:numcell-1);                    
                end
                
                while zone.trans(jj+1)>1.4 || trans.vec(end) > zone.dz(jj+1)
                numcell=numcell+1;
                zone.trans(jj+1)= get_factor(numcell,transtart,newdist);

                %Computes new transition factor
                trans.vec= transtart*zone.trans(jj+1).^(0:numcell-1);
                end
                
                while zone.trans(jj+1)<1.0 || zone.dz(jj+1)/trans.vec(end)>1.4
                numcell=numcell-1;
                zone.trans(jj+1)=get_factor(numcell,transtart,newdist);
                trans.vec= transtart * zone.trans(jj + 1).^(0:numcell-1);               
                end
                               
            end
            end

            zonemesh=[expan.vec(1:end) trans.vec(1:end)];
        
            
            mesh=[mesh zonemesh];
            trackzone=[trackzone ones(1,length(zonemesh))*jj];

            %Update next zone depth
            node(jj+1)=node(jj)-sum(zonemesh); 
            zone.depth(jj+1)=zone.depth(jj)+zone.depth(jj+1)-sum(zonemesh);
            zone.depth(jj)=sum(zonemesh); 
%             zone.dz(jj+1)=zonemesh(end);
            
        else
%             zone.dz(jj+1)=buffer(end);
            zonemesh=expan.vec(1:end);
            mesh=[mesh zonemesh];
            trackzone=[trackzone ones(1,length(zonemesh))*jj];

        end
        
% Back to general case 2  
        else

            %If the best overlapping cell is larger than the maxcellsize
            if expan.vec(flag.expan)>zone.sizemax(jj) || trans.vec(flag.trans)>zone.sizemax(jj)
            [zonemesh,zone,node]=maxcellsize(trans,expan,zone,flag,node,static(jj),jj);
            mesh=[mesh zonemesh];
            trackzone=[trackzone ones(1,length(zonemesh))*jj];
            jj=jj+1;
            continue
            else        

                zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
                trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);    

                [flag,expan,trans,zone,static,numcell,newdist]=fixtrans...
                    (flag,expan,trans,zone,static,numcell,newdist,jj);

                % If the transition is still too step, then the expan flag
                % is moved down of one until smooth transition.   
                while (zone.trans(jj+1)>1.4 || zone.trans(jj+1)<1.0 ||...
                        expan.vec(flag.expan)/trans.vec(end)>1.4) && flag.expan>2

                    flag.expan=flag.expan-1;
                    newdist=zone.depth(jj)-expan.length(flag.expan);
                    numcell=flag.trans+1;

                    zone.trans(jj+1)=get_factor(numcell,zone.dz(jj+1),newdist);
                    trans.vec=zone.dz(jj+1)*zone.trans(jj+1).^(0:numcell-1);

                    [flag,expan,trans,zone,static,numcell,newdist]=fixtrans...
                        (flag,expan,trans,zone,static,numcell,newdist,jj);
                end


            end

                zonemesh=[expan.vec(1:flag.expan) trans.vec(end:-1:1)];
                mesh=[mesh zonemesh];
                trackzone=[trackzone ones(1,length(zonemesh))*jj];
    %             sprintf('Zone %i starts at %6.3f',jj+1,zmax-sum(mesh))
        
        end
    end
end
    clear lookup
    jj=jj+1;
end
% end


%% Final Test - To see if the transitions are all <1.4. If not, the program 
% averages around the sharp corner. The direction of smoothing depends on
% the transition on both sides. Cells at the boundary of a zone are
% ignored. The script iterates until there is no sharp corners, or the
% corners are at the boundary.

mcell = length(mesh);
maxtrans = 1.4;

stride=mesh(1);
for ii=2:length(mesh)
    stride(ii)=stride(ii-1)+mesh(ii);
end
model=[zmax zmax-stride];
   

%Computes transitions everywhere
gradcell=mesh(1:end-1)./mesh(2:end);
gradcell(gradcell<=1)=1./gradcell(gradcell<=1);

corners = gradcell > maxtrans;
% Compute cells exclusion vector for boundary cells. This vector can be
% easily modified if we decide to relax the zone constraint.

exclusion = ones(1,mcell-1);
exclusion(1) = 0;

for ii=1:mcell-1
    if trackzone(ii) ~= trackzone(ii+1)
       exclusion(ii) = 0;
%        exclusion(ii+1) = 0;
    end
end
exclusion(end) = 0;

% while sum( corners .* exclusion ) > 0
%     for ii = find(corners .* exclusion == 1);
%         
%         if mesh(ii) > mesh(ii+1)
%            transup = mesh(ii+1) * maxtrans;
%            gap = mesh(ii) - transup;
%            mesh(ii+1) = mesh(ii+1) + gap/2;
%            mesh(ii) = transup + gap/2;
%         else
%            transup = mesh(ii) * maxtrans;
%            gap = mesh(ii+1) - transup;
%            mesh(ii) = mesh(ii) + gap/2;
%            mesh(ii+1) = transup + gap/2;
%         end
%     end
%     %Computes transitions everywhere
%     gradcell=mesh(1:end-1)./mesh(2:end);
%     gradcell(gradcell<=1)=1./gradcell(gradcell<=1);
% 
%     corners = gradcell > maxtrans;
% end
        
for ii=find( gradcell > maxtrans)
    sprintf('Found a transition > 1.4 at cell: %i in zone: %i', ii,trackzone(ii))
end


figure;imagesc(mesh);
sprintf('Final mesh is %6.2f m long',sum(mesh))
sprintf('Residual between VOI and mesh: %6.3f',node(1)-node(end)-sum(mesh))
write2UBC(zmax,mesh);
