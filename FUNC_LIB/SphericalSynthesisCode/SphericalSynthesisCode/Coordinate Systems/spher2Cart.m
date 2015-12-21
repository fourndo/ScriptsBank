function cartPoints=spher2Cart(points)
%%SPHER2CART  Convert points from spherical coordinates to Cartesian
%             coordinates
%
%INPUTS:  points  One or more points given in terms of range, azimuth and
%                 elevation, with the angles in radians. To convert
%                 N points, points is a 3XN matrix with each column
%                 having the format [range;azimuth; elevation].
%
%OUTPUTS:   cartPoints For N points, cartPoints is a 3XN matrix of the
%                      converted points with each column having the format
%                      [x;y;z].
%
%The conversion from spherical to Cartesian coordinates is given in
%R. L. Duncombe, "Computational techniques," in Explanatory Supplement
%to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
%Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
%ch. 14.4.4.1.
%
%Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
%is the angle above the x-y plane.
%
%December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

%Extract the coordinates
r=points(1,:);
azimuth=points(2,:);
elevation=points(3,:);

x=r.*cos(azimuth).*cos(elevation);
y=r.*sin(azimuth).*cos(elevation);
z=r.*sin(elevation);

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