/**
 * Receive an array of input parameters
 */
#include "../pvProj/ElectricalCalculation/MonoDescPWL.h"
#include "../pvProj/ModuleDefinition/OneDiodeModel.h"
#include "../pvProj/Utils/pvcommutils.h"
#include "../pvProj/testUtils.h"
#include "mex.h"
#include <omp.h>
#include <sched.h>
#include <string.h>

using namespace Eigen;
class MonoDescPWL;
class OneDiodeModel;

void mexVecODMpwlapproxAddseries(struct Params* inputParam, double* tol, double* lim,  double* w, unsigned int inputSize, unsigned int tolSize, unsigned int limSize, unsigned int wSize, double* outarray1, double* outarray2, unsigned int& length)
{
#ifdef TIMER
    double sXMLTime, eXMLTime;
    TIMESTAMP( sXMLTime );
#endif

    GetSimOptFile* filePtr = GetSimOptFile::Instance();

    char *RootPath = std::getenv("PVCPLUSPLUS_ROOT");
    if(RootPath == NULL)
    {
        std::cerr<<"The environment variable PVCPLUSPLUS_ROOT (the absolute path of folder pvCplusplus) should be preset!\n";
        std::cerr<<"Please read readme to get more detailed information about PVCPLUSPLUS_ROOT!\n";
        abort();
    }

    char SrcPath[PATH_MAX] = "\0"; // Absolute path of the resource file
    strcpy(SrcPath,RootPath);
    strcat(SrcPath,"/pvProj/Resources/SimOption.xml");

    filePtr->openOptionFile(SrcPath);
    filePtr->readOptionFile();

#ifdef TIMER
    TIMESTAMP( eXMLTime );
    std::cout<<"Pure XML read overhead...";
    std::cout<<"Elapsed time is "<< eXMLTime - sXMLTime <<" seconds\n";
#endif
    int i,j;
    std::vector<double> vecTol;
    std::vector<double> vecLim;
    std::vector<double> vecWeight;
      
    std::vector<MonoDescPWL> pwlset((inputSize));
    std::vector<MonoDescPWL*> set((inputSize));
    MonoDescPWL res;
    
    if (tolSize != 0)
        vecTol.assign(tol, tol+tolSize);
    if (limSize != 0)
        vecLim.assign(lim, lim+limSize);
    if (wSize != 0)
        vecWeight.assign(w, w+wSize);
         
/* vecLim[0] = lim[0];   
    vecLim[1] = lim[1];
    vecLim[3] = lim[2];
    vecLim[2] = lim[3]; */
  
#ifdef TIMER
    
    double sPWLApproxTime, ePWLApproxTime;
    TIMESTAMP( sPWLApproxTime );
#endif
    
    OneDiodeModel * output = new OneDiodeModel[inputSize];
    
    /// Vectorize.
#pragma omp parallel for default(shared) num_threads(4) if (inputSize > 6)
    for (i = 0; i < inputSize; i++)
    {
        output[i] = OneDiodeModel(inputParam[i]);
        output[i].pwlApprox(vecTol, vecLim);
        if (wSize != 0)
            { output[i].approxPointSet.col(0) *= vecWeight[i];  }
        pwlset[i] = MonoDescPWL(output[i].approxPointSet, {0.0});
        set[i] = &pwlset[i];
    }
   
#ifdef TIMER
    TIMESTAMP( ePWLApproxTime );
    std::cout<<"Pure C++ approx overhead...";
    std::cout<<"Elapsed time is "<< ePWLApproxTime - sPWLApproxTime <<" seconds\n";

    double sOutputTime, eOutputTime;
    TIMESTAMP( sOutputTime );
#endif
    /// Combined with addseries.
    res = MonoDescPWL(set, Dim::hori, vecTol, vecLim);
    length = res.pointSet.rows();

    for (j = 0; j < length; j++)
    {
        outarray1[j] = res.pointSet(j,0);
        outarray2[j] = res.pointSet(j,1);
    }
    
    delete [] output;
    
#ifdef TIMERR
    TIMESTAMP( eOutputTime );
    std::cout<<"Output overhead...";
    std::cout<<"Elapsed time is "<< eOutputTime - sOutputTime <<" seconds\n";
#endif
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
#ifdef TIMER
    double sWholeTime, eWholeTime;
    TIMESTAMP( sWholeTime );
#endif
    mwIndex idx = 0;
    int ifield, nfiields;
    double *outMatrix1;
    double *outMatrix2;
    double *tol, *lim, *weight;

    mwSize i, j, maxElem, sizeTol, sizeLim, sizeW, sizeInput;
#ifdef TIMER
    double sMallocTime, eMallocTime;
    TIMESTAMP( sMallocTime );
#endif

    struct Params* inputParam;
    sizeInput = mxGetN(prhs[0]);
    
    /// A vector of input parameters.
    inputParam = (struct Params*)mxMalloc( sizeInput * sizeof(struct Params) );
    
    maxElem = 4000;
    
    sizeTol = sizeLim = sizeW = 0;
    if (nrhs >= 2)
    {
        tol = mxGetPr(prhs[1]);
        sizeTol = mxGetN(prhs[1]);
    }
    if (nrhs >= 3)
    {
        lim = mxGetPr(prhs[2]);
        sizeLim = mxGetN(prhs[2]);
    }
    if (nrhs == 4)
    {
        weight = mxGetPr(prhs[3]);
        sizeW = mxGetM(prhs[3]);
    }
       
    for (i = 0; i < sizeInput; i++)
    {
        inputParam[i].Iph = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,0)));
        inputParam[i].Io = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,1)));
        inputParam[i].Rs = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,2)));
        inputParam[i].nVth = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,3)));
        inputParam[i].Rsh = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,4)));
        inputParam[i].di2Mutau = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,5)));
        inputParam[i].Vbi = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,6)));
        inputParam[i].Brev = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],i,7)));
    }
    
    /// Output.
    const mwSize dims = 1;
    plhs[0] = mxCreateDoubleMatrix(1, maxElem, mxREAL);
    outMatrix1 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, maxElem, mxREAL);
    outMatrix2 = mxGetPr(plhs[1]);
    plhs[2] = mxCreateNumericArray(1,&dims,mxUINT32_CLASS,mxREAL);
#ifdef TIMER
    TIMESTAMP( eMallocTime );
    std::cout<<"Malloc overhead...";
    std::cout<<"Elapsed time is "<< eMallocTime - sMallocTime <<" seconds\n"; 

   double sCallTime, eCallTime;
   TIMESTAMP( sCallTime );
#endif
    
    mexVecODMpwlapproxAddseries(inputParam, tol, lim, weight, sizeInput, sizeTol, sizeLim, sizeW, outMatrix1, outMatrix2, *(unsigned int*)mxGetData(plhs[2]));
#ifdef TIMER
    TIMESTAMP( eCallTime );
    std::cout<<"Call overhead...";
    std::cout<<"Elapsed time is "<< eCallTime - sCallTime <<" seconds\n"; 
#endif
    mxFree(inputParam);
#ifdef TIMER
    TIMESTAMP( eWholeTime );
    std::cout<<"Whole overhead...";
    std::cout<<"Elapsed time is "<< eWholeTime - sWholeTime <<" seconds\n";
#endif
}
