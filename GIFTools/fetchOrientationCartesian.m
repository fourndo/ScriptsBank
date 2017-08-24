% function ore = fetchOrientationCartesian(obj)
% Function converts angles [theta, alpha] to the closest 
% cartesian orientation. 
%
% Author: D.Fournier
% Date: 14 Aug 2017
%
% Modified from K. Davis
%

function ore = fetchOrientationCartesian(obj)
    % function ore = fetchOrientation(obj)
    % Returns 'x','y','z' orientation based on the flags
    % 1='x', 2='y', 3='z'
    theta = obj.getioData(:,{'Tx_ORIENTATION'});


    vert = theta < 45;
    rot = alpha > 45;

    if vert
        ore = 'z';
    else
        if rot
            ore = 'y';
        else
            ore = 'x';
        end
    end
end