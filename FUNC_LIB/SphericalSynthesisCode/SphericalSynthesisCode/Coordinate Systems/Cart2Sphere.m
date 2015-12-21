function points=Cart2Sphere(cartPoints)
%%CART2SPHERE Convert Cartesian coordinates to spherical coordinates.
%
%INPUTS: cartPoints A matrix of the points in Cartesian coordinates that
%                   are to be transformed into spherical coordinates. Each
%                   columns of cartPoints is of the format [x;y;z].
%
%OUTPUTS:   points  A matrix of the converted points. Each column of the
%                   matrix has the format [r;azimuth;elevation],
%                   with azimuth and elevation given in radians.
%
%The conversion from Cartesian to spherical coordinates is given in
%R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
%ch. 14.4.4.2.
%
%Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
%is the angle above the x-y plane.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Extract the coordinates
x=cartPoints(1,:);
y=cartPoints(2,:);
z=cartPoints(3,:);

r=sqrt(x.^2+y.^2+z.^2);
azimuth=atan2(y,x);
elevation=asin(z./r);

points=[r;azimuth;elevation];

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