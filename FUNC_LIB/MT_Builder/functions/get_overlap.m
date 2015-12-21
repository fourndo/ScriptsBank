function [flag_master,flag_slave]=get_overlap(master,slave,master_marker,slave_marker)
% Function to find the overlapping cell 
% [money]=overlap(master,slave)
% Input
% master: vector from which an index is found
% slave: vector being compared to
% master_marker:
% slave_marker:
%
% Output
% money: index of the best overlapping cell

        for ii=1:length(master_marker);
            dist=abs(master_marker(ii)-slave_marker);   %Distance to each cell
            nearest(ii)=find(dist==min(dist),1);          %Find closest cell      
            lookup(ii)=abs(master(ii)-slave(nearest(ii)));    %Compare cell sizes
        end
        index=find(lookup==min(lookup),1);
        flag_master=index;
        
        %Find the index of slave
        flag_slave=nearest(index);