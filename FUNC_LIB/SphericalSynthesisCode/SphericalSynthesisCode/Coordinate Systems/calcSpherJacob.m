function J=calcSpherJacob(point)
%%CALCSPHERJACOB  Compute the Jacobian matrix for a point in spherical
%                 [range;azimuth;elevation] coordinates.
%
%INPUTS: point   A point in the format [range;azimuth;elevation], where the
%                two angles are given in radians.
%
%OUTPUTS: J     The 3X3 Jacobian matrix. Each row is a components of range,
%               azimuth and elevation (in that order by row) with
%               derivatives taken with respect to [x,y,z] by column.
%
%The derivatives can be computed in a straightforward manner from
%the basic relation between spherical and Cartesian coordinates, which is
%given in
%R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
%ch. 14.4.4.1.
%among other sources.
%
%Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
%is the angle above the x-y plane. Note that singularities exist at the
%poles; that is when the elevation is +/-(pi/2).
%
%December 2013 David F.Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

CartPoint=spher2Cart(point);
x=CartPoint(1);
y=CartPoint(2);
z=CartPoint(3);

r=point(1);
J=zeros(3,3);

%Derivatives with respect to x.
J(1,1)=x/r;
J(2,1)=-y/(x^2+y^2);
J(3,1)=-x*z/(r^2*sqrt(x^2+y^2));

%Derivatives with respect to y.
J(1,2)=y/r;
J(2,2)=x/(x^2+y^2);
J(3,2)=-y*z/(r^2*sqrt(x^2+y^2));

%Derivatives with respect to z.
J(1,3)=z/r;
J(2,3)=0;
J(3,3)=sqrt(x^2+y^2)/r^2;

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