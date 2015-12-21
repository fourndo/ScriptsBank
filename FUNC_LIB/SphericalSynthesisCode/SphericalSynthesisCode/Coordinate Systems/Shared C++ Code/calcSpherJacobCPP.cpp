/*CALCSPHERJACOBCPP  A C++ function to compute the Jacobian matrix for
 *                   spherical coordinates.
 *
 *INPUTS: J    A pointer to an array of doubles with 9 elements to hold the
 *             result. The result is stored by column. The first column
 *             has derivatives with respect to x, the second y and the
 *             third z. The rows correspond to spherical radius, azimuth
 *             and elevation. 
 *     point   The 3X1 Cartesian point at which the Jacobian is to be
 *             evaluated.
 *OUTPUTS: None. The results are placed in J.
 *
 *Further comments are given in the Matlab file calcSpherJacob.
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
 **/

#include "CoordFuncs.hpp"
#include <math.h>

void calcSpherJacobCPP(double *J, const double *point) {
    double CartPoint[3],x,y,z,r;
    double x2y2Sum,xyDist,r2;
    
    spher2CartCPP(CartPoint,point);
    x=CartPoint[0];
    y=CartPoint[1];
    z=CartPoint[2];

    r=point[0];
    x2y2Sum=x*x+y*y;
    xyDist=sqrt(x2y2Sum);
    r2=r*r;
    
    //Derivatives with respect to x.
    J[0]=x/r;
    J[1]=-y/x2y2Sum;
    J[2]=-x*z/(r2*xyDist);

    //Derivatives with respect to y.
    J[3]=y/r;
    J[4]=x/x2y2Sum;
    J[5]=-y*z/(r2*xyDist);

    //Derivatives with respect to z.
    J[6]=z/r;
    J[7]=0;
    J[8]=xyDist/r2;
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
