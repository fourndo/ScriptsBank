/**SPHER2CART  A C++ implementation to convert points from spherical
*             coordinates to Cartesian coordinates
*
*INPUTS:  points  One or more points given in terms of range, azimuth and
*                 elevation, with the angles in radians. To convert
*                 N points, points is a 3XN matrix with each column
*                 having the format [range;azimuth; elevation].
*
*OUTPUTS:   cartPoints For N points, cartPoints is a 3XN matrix of the
*                      converted points with each column having the format
*                      [x;y;z].
*
*The conversion from spherical to Cartesian coordinates is given in
*R. L. Duncombe, "Computational techniques," in Explanatory Supplement
*to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
*Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
*ch. 14.4.4.1.
*
*Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
*is the angle above the x-y plane.
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*cartPoints=spher2Cart(points);
*
*January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
*/

#include <stddef.h>
//#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    size_t i, N;
    double *points,*retData;
    mxArray *retMat;
    
    if(nrhs!=1) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    if(mxGetM(prhs[0])!=3) {
       mexErrMsgTxt("The point has the wrong dimensionality.");
       return;
    }
    
    N=mxGetN(prhs[0]);//The number of points.
    checkRealDoubleArray(prhs[0]);
    points=(double*)mxGetData(prhs[0]);

    retMat=mxCreateDoubleMatrix(3, N,mxREAL);
    retData=(double*)mxGetData(retMat);
    
    for(i=0;i<N;i++) {
        spher2CartCPP(retData+3*i,points+3*i);
    }
    
    plhs[0]=retMat;
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
