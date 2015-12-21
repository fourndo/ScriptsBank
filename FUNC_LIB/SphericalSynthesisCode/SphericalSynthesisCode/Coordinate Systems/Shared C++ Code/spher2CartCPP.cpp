/*SPHER2CART  A C++ function to convert a point from 
 *
 *INPUTS: cartPoint    A pointer to an array of doubles with 3 elements to
 *                     hold the result in [x,y,z] order.
 *        point        The 3X1 point in spherical coordinates to be
 *                     converted, ordered [range;azimuth;elevation].
 *
 *OUTPUTS: None. The results are placed in cartPoint.
 *
 *Further comments are given in the Matlab file spher2Cart.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
 **/

//For sin and cos
#include <math.h>
#include "CoordFuncs.hpp"

void spher2CartCPP(double *cartPoint,const double *point) {
    double r, azimuth, elevation;
    double cosEl;
    r=point[0];
    azimuth=point[1];
    elevation=point[2];
    
    cosEl=cos(elevation);

    cartPoint[0]=r*cos(azimuth)*cosEl;
    cartPoint[1]=r*sin(azimuth)*cosEl;
    cartPoint[2]=r*sin(elevation);
}

/*LICENSE:
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
%OF RECIPIENT IN THE USE OF THE SOFTWARE.*/
