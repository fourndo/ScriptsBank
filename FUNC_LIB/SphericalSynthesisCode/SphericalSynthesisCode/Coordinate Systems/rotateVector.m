function vRot=rotateVector(v,u,theta)
%%ROTATEVECTOR  Given a three-dimensional vector v and a unit vector u,
%               rotate the vector an angle of theta counterclockwise about
%               the vector u.
%
%INPUTS:        v   The 3D vector that is to be rotated.
%               u   A unit vector representing the axis about which v is
%                   to be rotated.
%             theta The angle in radians by which v is to be rotated
%                   about u. The rotation angle is clockwise when one is
%                   looking in the same direction that the rotation axis
%                   points.
%
%OUTPUTS:       vRot The rotated vector.
%
%This simply implements the Rodrigues' rotation formula. The formula is
%given in equations (96) and (97) in
%M. D. Shuster, "A survey of attitude representations," The Journal of the
%Astronautical Sciences, vol. 41, no. 4, pp. 439-517, Oct.-Dec. 1993.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

vRot=v*cos(theta)+cross(u,v)*sin(theta)+u*dot(u,v)*(1-cos(theta));
end

%LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.