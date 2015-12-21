/*NALEGENDRECOSRAT A C++ implementation to evaluate
 *                 \bar{P}_{nm}(cos(theta))/u^m for all n from 0 to M and
 *                  for each n for all m from 0 to n. Also evaluate
 *                  D{\bar{P}_{nm}(cos(theta))}/u^m and
 *                  D2{\bar{P}_{nm}(cos(theta))}/u^m where D{} is the first
 *                  derivative operator with respect to theta and D2{} is
 *                  the second derivative operator with respect to theta.
 *                  \bar{P}_{nm}(x) is the fully normalized associated
 *                  Legendre function of x of degree n and order m.
 *                  u=sin(theta). All of the values can be scaled by a
 *                  factor of scalFac, if desired, to help prevent
 *                  overflows with high degrees and orders.
 *
 *INPUTS: theta     An angle in radians.
 *       maxDeg     The maximum degree and order of the output. This should
 *                  be at least 2.
 *     scalFactor   A scale factor to help prevent overflow of the results.
 *                  In the Holmes and Featherstone paper, discussed below,
 *                  a value of 10^(-280) is used.
 *
 *OUTPUTS: PBarUVals An instance of the ClusterSet class such that
 *                   PBarUVals(n+1,m+1)=scalFac*\bar{P}_{nm}(cos(theta))/u^m.
 *  dPBarUValsdTheta An instance of the ClusterSet class such that
 *                   dPBarUValsdTheta(n+1,m+1)=scalFac*D{\bar{P}_{nm}(cos(theta))}/u^m
 *d2PBarUValsdTheta2 An instance of the ClusterSet class such that
 *                   dPBarUValsdTheta(n+1,m+1)=scalFac*D2{\bar{P}_{nm}(cos(theta))}/u^m
 *
 *The modified forward row (MFR) algorithm of 
 *S. A. Holmes and W. E. Featherstone, "A unified approach to the Clenshaw
 *summation and the recursive computation of very high degree and
 *order normalised associated Legendre functions," Journal of Geodesy,
 *vol. 76, no. 5, pp. 279-299, May 2002.
 *is used to compute PBarUVals and dPBarUValsdTheta. For
 *d2PBarUValsdTheta2, the algorithm of 
 *S. A. Holmes and W. E. Featherstone, "Short note: Extending simplified
 *high-degree synthesis methods to second latitude derivatives of
 *geopotential," Journal of Geodesy, vol. 76, no. 8, pp. 447-450, Nov.
 *2002.
 *is used. However, the paper contains a typo. The first term of the first
 *unnumbered equation should be multiplied by an additional (1/u). This
 *function uses the correct formulation.
 *
 *Additional comments are given in the Matlab implementation. This
 *implementation tends to be approximately 11,000 times faster than the
 *Matlab implementation.
 *
 *Note that this function requires a large amount of memory, since Matlab's
 *mxSetProperty function copies the data rather than allowing pointers to
 *be passed. Thus, more than double the normal amount of space is required
 *for the return values.
 *
 **The algorithm can be compiled for use in Matlab using the 
 *CompileCLibraries function.
 *
 *The algorithm is run in Matlab using the command format
 *[PBarUVals,dPBarUValsdTheta,d2PBarUValsdTheta2]=NALegendreCosRat(theta,M,scalFactor);
 *
 *January 2014 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
 */

#include "mex.h"
/* This header validates inputs and includes a header needed to handle
 * Matlab matrices.*/
#include "MexValidation.h"
#include "mathFuncs.hpp"

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[]) {
    double theta, scalFactor;
    ClusterSetCPP<double> PBarUVals;
    ClusterSetCPP<double> dPBarUValsdTheta;//The first derivatives
    ClusterSetCPP<double> d2PBarUValsdTheta2;//The second derivatives
    size_t M, numPBarU,i;
    mxArray *CSRetVal;
    mxArray *clusterElsMATLAB,*clusterSizesMATLAB, *offsetArrayMATLAB;
    
    if(nrhs!=3){
        mexErrMsgTxt("Incorrect number of inputs.");
        return;
    }

    theta=getDoubleFromMatlab(prhs[0]);
    M=getSizeTFromMatlab(prhs[1]);
    scalFactor=getDoubleFromMatlab(prhs[2]);
    
    if(M<2) {
       mexErrMsgTxt("The maximum order should be at least 2.");
       return; 
    }
    
    numPBarU=(M+1)*(M+2)/2;
    
    //Allocate space for the results.
    clusterElsMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
    clusterSizesMATLAB=allocUnsignedSizeMatInMatlab(M+1,1);
    offsetArrayMATLAB=allocUnsignedSizeMatInMatlab(M+1,1);
    
    PBarUVals.numClust=M+1;
    PBarUVals.totalNumEl=numPBarU;
    PBarUVals.clusterEls=(double*)mxGetData(clusterElsMATLAB);
    PBarUVals.offsetArray=(size_t*)mxGetData(offsetArrayMATLAB);
    PBarUVals.clusterSizes=(size_t*)mxGetData(clusterSizesMATLAB);
    
    //Initialize the offset array and cluster sizes.
    PBarUVals.offsetArray[0]=0;
    PBarUVals.clusterSizes[0]=1;
    for(i=1;i<=M;i++){
        PBarUVals.clusterSizes[i]=i+1;
        PBarUVals.offsetArray[i]=PBarUVals.offsetArray[i-1]+PBarUVals.clusterSizes[i-1];
    }
    
    NALegendreCosRatCPP(PBarUVals, theta,scalFactor);
    
    //Set the first return value
    mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
    mxSetProperty(CSRetVal,0,"clusterEls",clusterElsMATLAB);
    mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
    mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);
    
    plhs[0]=CSRetVal;

    if(nlhs>1) {//Compute the first derivatives, if they are desired.
        mxArray *clusterEls1stDerivMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
        
        dPBarUValsdTheta.numClust=M+1;
        dPBarUValsdTheta.totalNumEl=numPBarU;
        dPBarUValsdTheta.clusterEls=(double*)mxGetData(clusterEls1stDerivMATLAB);
        dPBarUValsdTheta.offsetArray=(size_t*)mxGetData(offsetArrayMATLAB);
        dPBarUValsdTheta.clusterSizes=(size_t*)mxGetData(clusterSizesMATLAB);

        NALegendreCosRatDerivCPP(dPBarUValsdTheta, PBarUVals, theta);
        
        //Set the second return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls1stDerivMATLAB);
        mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
        mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);

        plhs[1]=CSRetVal;
        mxDestroyArray(clusterEls1stDerivMATLAB);
    }
    
    if(nlhs>2) {//Compute the second derivatives if they are desired.
        mxArray *clusterEls2ndDerivMATLAB=mxCreateDoubleMatrix(numPBarU,1,mxREAL);
        
        d2PBarUValsdTheta2.numClust=M+1;
        d2PBarUValsdTheta2.totalNumEl=numPBarU;
        d2PBarUValsdTheta2.clusterEls=(double*)mxGetData(clusterEls2ndDerivMATLAB);
        d2PBarUValsdTheta2.offsetArray=(size_t*)mxGetData(offsetArrayMATLAB);
        d2PBarUValsdTheta2.clusterSizes=(size_t*)mxGetData(clusterSizesMATLAB);
        
        NALegendreCosRatDeriv2CPP(d2PBarUValsdTheta2, dPBarUValsdTheta, PBarUVals,theta);
        
        //Set the third return value
        mexCallMATLAB(1,&CSRetVal,0, 0, "ClusterSet");
        mxSetProperty(CSRetVal,0,"clusterEls",clusterEls2ndDerivMATLAB);
        mxSetProperty(CSRetVal,0,"clusterSizes",clusterSizesMATLAB);
        mxSetProperty(CSRetVal,0,"offsetArray",offsetArrayMATLAB);

        plhs[2]=CSRetVal;
        mxDestroyArray(clusterEls2ndDerivMATLAB);
    }
    
    //Free the buffers. The mxSetProperty command copied the data.
    mxDestroyArray(clusterElsMATLAB);
    mxDestroyArray(clusterSizesMATLAB);
    mxDestroyArray(offsetArrayMATLAB);
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
