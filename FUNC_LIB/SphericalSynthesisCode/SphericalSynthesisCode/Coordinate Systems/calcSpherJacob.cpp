/*CALCSPHERJACOB  A C++ function to compute the Jacobian matrix for a point
*                 in spherical [range;azimuth;elevation] coordinates.
*
*INPUTS: point   A point in the format [range.azimuth;elevation], where the
*                two angles are given in radians.
*
*OUTPUTS: J     The 3X3 Jacobian matrix. Each row is a component of range,
*               azimuth and elevation (in that order by row) with
*               derivatives taken with respect to [x,y,z] by column.
*
*The derivatives can be computed in a straightforward manner from
*the basic relation between spherical and Cartesian coordinates, which is
*given in
*R. L. Duncombe, "Computational techniques," in Explanatory Supplement
*to the Astronomical Almanac, 3rd ed., S. E. Urban and P. K.
*Seidelmann, Eds. Mill Valley, CA: University Science Books, 2013,
*ch. 14.4.4.1.
*among other sources.
*
*Azimuth is an angle measured from the x-axis in the x-y plane. Elevation
*is the angle above the x-y plane. Note that singularities exist at the
*poles; that is when the elevation is +/-(pi/2).
*
*The algorithm can be compiled for use in Matlab  using the 
*CompileCLibraries function.
*
*The algorithm is run in Matlab using the command format
*J=calcSpherJacob(point);
*
*January 2014 David F.Crouse, Naval Research Laboratory, Washington D.C.
*(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
*/

#include "matrix.h"
#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "CoordFuncs.hpp"
        
void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double *point, *retData;
    mxArray *retMat;
    
    if(nrhs!=1) {
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }
    
    if(nlhs>1) {
        mexErrMsgTxt("Wrong number of outputs.");
        return;
    }
    
    
    if(mxGetM(prhs[0])!=3||mxGetN(prhs[0])!=1) {
       mexErrMsgTxt("The point has the wrong dimensionality.");
       return;
    }
    
    checkRealDoubleArray(prhs[0]);
    point=(double*)mxGetData(prhs[0]);

    retMat=mxCreateDoubleMatrix(3, 3,mxREAL);
    retData=(double*)mxGetData(retMat);
    
    calcSpherJacobCPP(retData, point);
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
