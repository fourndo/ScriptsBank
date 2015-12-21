function cartPoints=ellips2Cart(points,a,f)
%%ELLIPS2CART Convert ellipsoidal coordinates to ECEF Cartesian
%             coordinates.
%
%INPUTS:    points  One or more points given in geodetic latitude and
%                   longitude, in radians, and height, in meters that are
%                   to be converted to Cartesian coordinates. To convert
%                   N points, points is a 3XN matrix with each column
%                   having the format [latitude;longitude; height].
%           a       The semi-major axis of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84SemiMajorAxis is used.
%           f       The flattening factor of the reference ellipsoid. If
%                   this argument is omitted, the value in
%                   Constants.WGS84Flattening is used.
%
%OUTPUTS:   cartPoints For N points, cartPoints is a 3XN matrix of the
%                      converted points with each column having the format
%                      [x;y;z].
%
%The conversions are mentioned in the paper "Simulating Aerial Targets in
%3D Accounting for the Earth's Curvature" by David F. Crouse.
%
%September 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

if(nargin<3)
    f=Constants.WGS84Flattening;
end

if(nargin<2)
    a=Constants.WGS84SemiMajorAxis;
end

%The geodetic latitudes
phi=points(1,:);
%The longitudes
lambda=points(2,:);
%The altitudes
h=points(3,:);

sinP=sin(phi);
cosP=cos(phi);
sinL=sin(lambda);
cosL=cos(lambda);

%The square of the first numerical eccentricity
e2=2*f-f^2;
%The normal radii of curvature.
Ne=a./sqrt(1-e2*sinP.^2);

x=(Ne+h).*cosP.*cosL;
y=(Ne+h).*cosP.*sinL;
z=(Ne*(1-e2)+h).*sinP;

cartPoints=[x;y;z];
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