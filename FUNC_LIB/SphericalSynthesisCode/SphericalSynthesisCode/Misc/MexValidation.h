/* This is a suite of C and C++ functions for type checking inputs to mex
 * files and performing other helper actions when interfacing Matlab and C
 * or C++ code.
 *
 *The following functions that can be used in both C and C++ programs:
 *checkRealDoubleArray      Gives an error if an array is not composed of
 *                          real doubles.
 *verifySizeReal            Verify that a Matlab matrix is 2D, has the
 *                          desired dimensions, and is real.
 *convert2DReal2DoubleMat   Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          doubles. The returned Matrix should be freed
 *                          with mxDestroyArray. the input matrix is not
 *                          modified.
 *convert2DReal2SignedIntMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          ints. The returned Matrix should be freed
 *                          with mxDestroyArray. the input matrix is not
 *                          modified.
 *convert2DReal2UnsignedIntMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          unsigned ints. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *convert2DReal2SignedSizeMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          ptrdiff_t. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *convert2DReal2UnsignedSizeMat Convert a Matlab matrix of real values for
 *                          various Matlab data types to a Matlab matrix of
 *                          size_t. The returned Matrix should be
 *                          freed with mxDestroyArray. the input matrix is
 *                          not modified.
 *getIntFromMatlab          Turns a passed Matlab data type into a real
 *                          integer. Provides an error if the Matlab data
 *                          is of a type that can not be transformed into
 *                          an integer or if it is complex.
 *getSizeTFromMatlab        Turns a passed Matlab data type into a
 *                          positive, real size_t integer value. Provides
 *                          an error if the Matlab data is of a type that
 *                          can not be transformed into an size_t value or
 *                          if it is negative or complex.
 *copySizeTArrayFromMatlab  Returns a copy of an array from Matlab as an
 *                          array of size_t values and indicated the length
 *                          of the array by modifying an input parameter.
 *getDoubleFromMatlab       Turns a passed Matlab data type into a real
 *                          double variable. Provides an error if the
 *                          Matlab data is of a type that can not be
 *                          transformed into a double or if it is
 *                          complex.
 *getBoolFromMatlab         Turns a passed Matlab data type into a real
 *                          boolean variable. Provides an error if the
 *                          Matlab data is of a type that can not be
 *                          transformed into a boolean or if it is
 *                          complex.
 *doubleMat2Matlab          Convert an array or matrix of doubles into
 *                          a matrix to be returned to Matlab.
 *allocUnsignedSizeMatInMatlab Allocate a matrix for holding size_t values
 *                          that can be used in Matlab.
 *allocSignedSizeMatInMatlab Allocate a matrix for holding ptrdiff_t values
 *                          that can be used in Matlab.
 *allocUnsignedIntMatInMatlab Allocate a matrix for holding unsigned int
 *                          values that can be used in Matlab.
 *allocSignedIntMatInMatlab Allocate a matrix for holding int values
 *                          that can be used in Matlab.
 *unsignedSizeMat2Matlab    Convert an array or matrix of size_t values
 *                          into a matrix to be returned by Matlab.
 *signedSizeMat2Matlab      Convert an array or matrix of ptrdiff_t values
 *                          into a matrix to be returned to Matlab.
 *getScalarMatlabClassConst Get a double scalar constant that is defined in
 *                          a Matlab class.
 *
 *The following functions are for use when interfacing C and C++ programs:
 *ptr2Matlab                Convert a pointer to a class instance into a
 *                          data type that can be returned to Matlab so
 *                          that the class can be used during another call
 *                          to a mex function.
 *Matlab2Ptr                Convert a Matlab value obtained from ptr2Matlab
 *                          back to a pointer to a C++ class.
 *
 *DEPENDENCIES: matrix.h
 *              mex.h
 *              algorithms if using C++
 *              cstfdef if using C++
 *              string.h if using C
 *              stddef.h if using C
 *              stdbool.h if using C99 or later.
 *
 *December 2013 David F. Crouse, Naval Research Laboratory, Washington D.C.
 *(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.
 */

#ifndef MEXHELP
#define MEXHELP
#include "matrix.h"
#include "mex.h"

#ifdef __cplusplus
//For memcpy
#include <algorithm>
//Defines the size_t and ptrdiff_t types
#include <cstddef>
#else
//For memcpy
#include <string.h>
//Defines the size_t and ptrdiff_t types
#include <stddef.h>

//This is needed for the bool type to be defined in C.
//C99 has stdbool, earlier versions do not.
#if __STDC_VERSION__>=199901L
#include <stdbool.h>
#else
//This is for the bool that Matlab's headers might define.
#ifndef _bool_T
#define false 0
#define true 1
#define bool int
#endif
#endif
#endif

void checkRealDoubleArray(const mxArray *val);
void verifySizeReal(size_t M, size_t N, const mxArray *val);
mxArray *convert2DReal2DoubleMat(const mxArray *val);
mxArray *convert2DReal2SignedIntMat(const mxArray *val);
mxArray *convert2DReal2UnsignedIntMat(const mxArray *val);
mxArray *convert2DReal2SignedSizeMat(const mxArray *val);
mxArray *convert2DReal2UnsignedSizeMat(const mxArray *val);

int getIntFromMatlab(const mxArray *val);
size_t getSizeTFromMatlab(const mxArray *val);
size_t *copySizeTArrayFromMatlab(const mxArray *val, size_t *arrayLen);
double getDoubleFromMatlab(const mxArray *val);
bool getBoolFromMatlab(const mxArray *val);

mxArray *doubleMat2Matlab(const double *arr,const size_t numRow, const size_t numCol);
mxArray *allocUnsignedSizeMatInMatlab(const size_t numRow, const size_t numCol);
mxArray *allocSignedSizeMatInMatlab(const size_t numRow, const size_t numCol);
mxArray *allocUnsignedIntMatInMatlab(const size_t numRow, const size_t numCol);
mxArray *allocSignedIntMatInMatlab(const size_t numRow, const size_t numCol);
mxArray *unsignedSizeMat2Matlab(const size_t *arr, const size_t numRow, const size_t numCol);
mxArray *signedSizeMat2Matlab(const ptrdiff_t *arr, const size_t numRow, const size_t numCol);
double getScalarMatlabClassConst(const char *className, const char *constName);

void checkRealDoubleArray(const mxArray *val){
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    if(mxGetClassID(val)!=mxDOUBLE_CLASS) {
        mexErrMsgTxt("A parameter that should be a real double is of a different data type.");
    }
}

void verifySizeReal(size_t M, size_t N, const mxArray *val) {
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix has complex components.");
    }
    
    if(mxGetM(val)!=M||mxGetN(val)!=N) {
        mexErrMsgTxt("A parameter has the wrong dimensionality.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
}

mxArray *convert2DReal2DoubleMat(const mxArray *val) {
    mxArray *retMat;
    double *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    retMat=mxCreateDoubleMatrix(M,N,mxREAL);
    retData=(double*)mxGetData(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *data=(char*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *data=(bool*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=(double*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=(float*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=(int8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=(uint8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=(int16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=(uint16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=(int32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=(uint32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=(int64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=(uint64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(double)data[i];
            }
            break;
        }
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *convert2DReal2SignedIntMat(const mxArray *val) {
    mxArray *retMat;
    int *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocSignedIntMatInMatlab(M,N);
    retData=(int*)mxGetData(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *data=(char*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *data=(bool*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=(double*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=(float*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=(int8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=(uint8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=(int16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=(uint16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=(int32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=(uint32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=(int64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=(uint64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(int)data[i];
            }
            break;
        }
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *convert2DReal2UnsignedIntMat(const mxArray *val) {
    mxArray *retMat;
    unsigned int *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocUnsignedIntMatInMatlab(M,N);
    retData=(unsigned int*)mxGetData(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *data=(char*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *data=(bool*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=(double*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=(float*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=(int8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=(uint8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=(int16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=(uint16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=(int32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=(uint32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=(int64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=(uint64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(unsigned int)data[i];
            }
            break;
        }
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *convert2DReal2SignedSizeMat(const mxArray *val) {
    mxArray *retMat;
    ptrdiff_t *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocSignedSizeMatInMatlab(M,N);
    retData=(ptrdiff_t*)mxGetData(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *data=(char*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *data=(bool*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=(double*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=(float*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=(int8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=(uint8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=(int16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=(uint16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=(int32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=(uint32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=(int64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=(uint64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(ptrdiff_t)data[i];
            }
            break;
        }
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

mxArray *convert2DReal2UnsignedSizeMat(const mxArray *val) {
    mxArray *retMat;
    size_t *retData;
    size_t M,N,numElements,i;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles has complex components.");
    }
    
    if(mxIsEmpty(val)) {
        mexErrMsgTxt("A parameter that should be real matrix of doubles is empty.");
    }
    
    if(mxGetNumberOfDimensions(val)>2){
        mexErrMsgTxt("A parameter has too many dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    numElements=M*N;
    
    retMat=allocUnsignedSizeMatInMatlab(M,N);
    retData=(size_t*)mxGetData(retMat);
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *data=(char*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *data=(bool*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *data=(double*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *data=(float*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *data=(int8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *data=(uint8_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *data=(int16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *data=(uint16_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *data=(int32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *data=(uint32_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *data=(int64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *data=(uint64_T*)mxGetData(val);
            for(i=0;i<numElements;i++) {
                retData[i]=(size_t)data[i];
            }
            break;
        }
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retMat;
}

int getIntFromMatlab(const mxArray *val) {
    int retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(int)*(char*)mxGetData(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(int)*(bool*)mxGetData(val);
            break;
        case mxDOUBLE_CLASS:
            retVal=(int)*(double*)mxGetData(val);
            break;
        case mxSINGLE_CLASS:
            retVal=(int)*(float*)mxGetData(val);
            break;
        case mxINT8_CLASS:
            retVal=(int)*(int8_T*)mxGetData(val);
            break;
        case mxUINT8_CLASS:
            retVal=(int)*(uint8_T*)mxGetData(val);
            break;
        case mxINT16_CLASS:
            retVal=(int)*(int16_T*)mxGetData(val);
            break;
        case mxUINT16_CLASS:
            retVal=(int)*(uint16_T*)mxGetData(val);
            break;
        case mxINT32_CLASS:
            retVal=(int)*(int32_T*)mxGetData(val);
            break;
        case mxUINT32_CLASS:
            retVal=(int)*(uint32_T*)mxGetData(val);
            break;
        case mxINT64_CLASS:
            retVal=(int)*(int64_T*)mxGetData(val);
            break;
        case mxUINT64_CLASS:
            retVal=(int)*(uint64_T*)mxGetData(val);
            break;
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

size_t getSizeTFromMatlab(const mxArray *val) {
    size_t retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(size_t)*(char*)mxGetData(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(size_t)*(bool*)mxGetData(val);
            break;
        case mxDOUBLE_CLASS:
        {
            double temp=*(double*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxSINGLE_CLASS:
        {
            float temp=*(float*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T temp=*(int8_T*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT8_CLASS:
            retVal=(size_t)*(uint8_T*)mxGetData(val);
            break;
        case mxINT16_CLASS:
        {
            int16_T temp=*(int16_T*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT16_CLASS:
            retVal=(size_t)*(uint16_T*)mxGetData(val);
            break;
        case mxINT32_CLASS:
        {
            int32_T temp=*(int32_T*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT32_CLASS:
            retVal=(size_t)*(uint32_T*)mxGetData(val);
            break;
        case mxINT64_CLASS:
        {
            int64_T temp=*(int64_T*)mxGetData(val);
            if(temp<0) {
                mexErrMsgTxt("A parameter that should be positive is not.");
            }
            retVal=(size_t)temp;
            break;
        }
        case mxUINT64_CLASS:
            retVal=(size_t)*(uint64_T*)mxGetData(val);
            break;
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

size_t *copySizeTArrayFromMatlab(const mxArray *val, size_t *arrayLen) {
//The length of the returned array is placed into the arrayLen variable.
    size_t M,N,i;
    size_t numEl=0;
    size_t *retVal=NULL;
    
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfDimensions(val)>2) {
        mexErrMsgTxt("A parameter that should be an array has extra dimensions.");
    }
    
    M=mxGetM(val);
    N=mxGetN(val);
    
    if(M==1&&N>=1) {
        numEl=N;
    } else if(M>1&&N==1){
        numEl=M;
    } else {
        mexErrMsgTxt("A parameter that should be an array has extra dimensions.");
    }
    
    *arrayLen=numEl;
    
    #ifdef __cplusplus
    retVal= new size_t[numEl];
    #else
    retVal=mxMalloc(numEl*sizeof(size_t));
    #endif     
            
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
        {
            char *mexData=(char*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxLOGICAL_CLASS:
        {
            bool *mexData=(bool*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxDOUBLE_CLASS:
        {
            double *mexData=(double*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxSINGLE_CLASS:
        {
            float *mexData=(float*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT8_CLASS:
        {
            int8_T *mexData=(int8_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT8_CLASS:
        {
            uint8_T *mexData=(uint8_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT16_CLASS:
        {
            int16_T *mexData=(int16_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT16_CLASS:
        {
            uint16_T *mexData=(uint16_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxINT32_CLASS:
        {
            int32_T *mexData=(int32_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT32_CLASS:
        {
            uint32_T *mexData=(uint32_T*)mxGetData(val);
            
            if(sizeof(size_t)==sizeof(uint64_T)) {
                memcpy(retVal,mexData,numEl*sizeof(uint32_T));
            } else {
                for(i=0;i<numEl;i++) {
                    retVal[i]=(size_t)mexData[i];
                }
            }
            break;
        }
        case mxINT64_CLASS:
        {
            int64_T *mexData=(int64_T*)mxGetData(val);
            
            for(i=0;i<numEl;i++) {
                if(mexData[i]<0) {
                    #ifdef __cplusplus
                    delete[] retVal;
                    #else
                    mxFree(retVal);
                    #endif
                    
                    mexErrMsgTxt("A parameter that should be positive is not.");
                }
                
                retVal[i]=(size_t)mexData[i];
            }
            break;
        }
        case mxUINT64_CLASS:
        {
            uint64_T *mexData=(uint64_T*)mxGetData(val);
            
            if(sizeof(size_t)==sizeof(uint64_T)) {
                memcpy(retVal,mexData,numEl*sizeof(uint64_T));
            } else {
                for(i=0;i<numEl;i++) {
                    retVal[i]=(size_t)mexData[i];
                }
            }
            break;
        }
        default:
            #ifdef __cplusplus
            delete[] retVal;
            #else
            mxFree(retVal);
            #endif  
            
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    
    return retVal;
}

double getDoubleFromMatlab(const mxArray *val) {
    double retVal=0;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=(double)*(char*)mxGetData(val);
            break;
        case mxLOGICAL_CLASS:
            retVal=(double)*(bool*)mxGetData(val);
            break;
        case mxDOUBLE_CLASS:
            retVal=(double)*(double*)mxGetData(val);
            break;
        case mxSINGLE_CLASS:
            retVal=(double)*(float*)mxGetData(val);
            break;
        case mxINT8_CLASS:
            retVal=(double)*(int8_T*)mxGetData(val);
            break;
        case mxUINT8_CLASS:
            retVal=(double)*(uint8_T*)mxGetData(val);
            break;
        case mxINT16_CLASS:
            retVal=(double)*(int16_T*)mxGetData(val);
            break;
        case mxUINT16_CLASS:
            retVal=(double)*(uint16_T*)mxGetData(val);
            break;
        case mxINT32_CLASS:
            retVal=(double)*(int32_T*)mxGetData(val);
            break;
        case mxUINT32_CLASS:
            retVal=(double)*(uint32_T*)mxGetData(val);
            break;
        case mxINT64_CLASS:
            retVal=(double)*(int64_T*)mxGetData(val);
            break;
        case mxUINT64_CLASS:
            retVal=(double)*(uint64_T*)mxGetData(val);
            break;
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

bool getBoolFromMatlab(const mxArray *val) {
    bool retVal=false;
    if(mxIsComplex(val)==true) {
        mexErrMsgTxt("A parameter that should be real has complex components.");
    }
    
    if(mxGetNumberOfElements(val)!=1) {
        mexErrMsgTxt("A parameter that should be scalar is not.");
    }
    
    switch(mxGetClassID(val)){
        case mxCHAR_CLASS:
            retVal=*(char*)mxGetData(val)!=0.0;
            break;
        case mxLOGICAL_CLASS:
            retVal=*(bool*)mxGetData(val)!=0;
            break;
        case mxDOUBLE_CLASS:
            retVal=*(double*)mxGetData(val)!=0.0;
            break;
        case mxSINGLE_CLASS:
            retVal=*(float*)mxGetData(val)!=0;
            break;
        case mxINT8_CLASS:
            retVal=*(int8_T*)mxGetData(val)!=0;
            break;
        case mxUINT8_CLASS:
            retVal=*(uint8_T*)mxGetData(val)!=0;
            break;
        case mxINT16_CLASS:
            retVal=*(int16_T*)mxGetData(val)!=0;
            break;
        case mxUINT16_CLASS:
            retVal=*(uint16_T*)mxGetData(val)!=0;
            break;
        case mxINT32_CLASS:
            retVal=*(int32_T*)mxGetData(val)!=0;
            break;
        case mxUINT32_CLASS:
            retVal=*(uint32_T*)mxGetData(val)!=0;
            break;
        case mxINT64_CLASS:
            retVal=*(int64_T*)mxGetData(val)!=0;
            break;
        case mxUINT64_CLASS:
            retVal=*(uint64_T*)mxGetData(val)!=0;
            break;
        default:
            mexErrMsgTxt("A parameter is of a data type that can not be used.");
    }
    return retVal;
}

mxArray *doubleMat2Matlab(const double *arr,const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    void *dataPtr;
    
    retMat=mxCreateDoubleMatrix(numRow,numCol,mxREAL);
    dataPtr=mxGetData(retMat);
    memcpy(dataPtr,arr,numRow*numCol*sizeof(double));

    return retMat;
}

mxArray *allocUnsignedSizeMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(size_t)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *allocSignedSizeMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(ptrdiff_t)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *allocUnsignedIntMatInMatlab(const size_t numRow, const size_t numCol) {
    mxArray *retMat=NULL;
    
    switch(sizeof(unsigned int)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxUINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *allocSignedIntMatInMatlab(const size_t numRow, const size_t numCol){
    mxArray *retMat=NULL;
    
    switch(sizeof(int)) {
        case 4:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retMat = mxCreateNumericMatrix(numRow,numCol,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    return retMat;
}

mxArray *unsignedSizeMat2Matlab(const size_t *arr, const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    void *dataPtr;
    
    retMat=allocUnsignedSizeMatInMatlab(numRow, numCol);
    
    dataPtr=mxGetData(retMat);
    memcpy(dataPtr,arr,numRow*numCol*sizeof(size_t));

    return retMat;
}

mxArray *signedSizeMat2Matlab(const ptrdiff_t *arr, const size_t numRow, const size_t numCol) {
    mxArray *retMat;
    void *dataPtr;
    
    retMat=allocSignedSizeMatInMatlab(numRow,numCol);
    
    dataPtr=mxGetData(retMat);
    memcpy(dataPtr,arr,numRow*numCol*sizeof(ptrdiff_t));

    return retMat;
}

double getScalarMatlabClassConst(const char *className, const char *constName) {
/*Whereas in Matlab one can access constants in classes without
 *instantiating the class, in mex, the class must be instantiated first.
 *Class name and const name are null-terminated strings. For example, if one wishes to access
 *Constants.standardTemp, then one can use
 *const char className[]="Constants";
 *const char constName[]="standardTemp";
 **/
    mxArray *constantClass, *valMATLAB;
    double retVal;

    mexCallMATLAB(1,&constantClass,0,NULL,className);//Load the Constants class.
    valMATLAB=mxGetProperty(constantClass, 0,constName);

    if(valMATLAB==NULL) {
        mexErrMsgTxt("An eror occurred while trying to access a class constant.");
    }
    mxDestroyArray(constantClass);
    retVal=getDoubleFromMatlab(valMATLAB);
    mxDestroyArray(valMATLAB);
    return retVal;
}

#ifdef __cplusplus

template<typename T>
mxArray *ptr2Matlab(T thePointer) {
/*Ptr2Matlab Convert a C++ pointer to a Matlab matrix so that it can be
 *           returned to Matlab between calls to mex functions.
 */
    mxArray *retArray=NULL;
    
    switch(sizeof(void*)) {
        case 4:
            retArray = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
            break;
        case 8:
            retArray = mxCreateNumericMatrix(1,1,mxINT64_CLASS,mxREAL);
            break;
        default:
            mexErrMsgTxt("The integer size of this computer is neither 64 nor 32 bit. Thus, the Matlab data type for pointer conversion could not be determined.");
    }
    
    *(void**)mxGetData(retArray)=reinterpret_cast<void*>(thePointer);
    return retArray;
}

template<typename T>
T Matlab2Ptr(const mxArray *matlabPtr){
/**MATLAB2PTR Convert a 1X1 Matlab matrix in which a pointer has been saved
 *            using ptr2Matlab back into a pointer for use in a mex file.
 */
    T thePtr;
    thePtr=reinterpret_cast<T>(*(void**)mxGetData(matlabPtr));    
    return thePtr;
}

#endif
#endif

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
