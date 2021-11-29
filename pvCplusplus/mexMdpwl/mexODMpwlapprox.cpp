/**
 * Receive an input parameter.
 */
#include "../pvProj/ModuleDefinition/OneDiodeModel.h"
#include "../pvProj/Utils/pvcommutils.h"
#include "../pvProj/testUtils.h"
#include "mex.h"
#include <omp.h>
#include <sched.h>
#include <string.h>

using namespace Eigen;
//namespace fs = boost::filesystem;
class OneDiodeModel;

void mexODMpwlapprox(struct Params inputParam, double* tol, double* lim, unsigned int tolSize, unsigned int limSize, double* outarray1, double* outarray2, unsigned int& length)
{
#ifdef TIMER
    double sXMLTime, eXMLTime;
    TIMESTAMP( sXMLTime );
#endif

    GetSimOptFile* filePtr = GetSimOptFile::Instance();

    char SrcPath[PATH_MAX] = "\0"; // Absolute path of the resource file
    realpath("./", SrcPath);
    
    char *substr = strstr(SrcPath, "solar-simulation");
    strcpy(substr, "solar-simulation/pvCplusplus/pvProj/Resources/SimOption.xml");

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
    
    if (tolSize != 0)
        vecTol.assign(tol, tol+tolSize);
    if (limSize != 0)
        vecLim.assign(lim, lim+limSize);
    
/* vecLim[0] = lim[0];   
    vecLim[1] = lim[1];
    vecLim[3] = lim[2];
    vecLim[2] = lim[3];*/
   
#ifdef TIMER
    double sPWLApproxTime, ePWLApproxTime;
    TIMESTAMP( sPWLApproxTime );
#endif
    
    /// Corresponds to onediodepwlapprox
    OneDiodeModel output(inputParam);
    output.pwlApprox(vecTol, vecLim);
#ifdef TIMER
    TIMESTAMP( ePWLApproxTime );
    std::cout<<"Pure C++ approx overhead...";
    std::cout<<"Elapsed time is "<< ePWLApproxTime - sPWLApproxTime <<" seconds\n";

    double sOutputTime, eOutputTime;
    TIMESTAMP( sOutputTime );
#endif
    length = output.approxPointSet.rows();
    
    /// Return the array to mex.
    for (j = 0; j < length; j++)
    {
        outarray1[j] = output.approxPointSet(j,0);
        outarray2[j] = output.approxPointSet(j,1);
    }
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
    double *tol, *lim;

    mwSize i, j, maxElem, sizeTol, sizeLim;
#ifdef TIMER
    double sMallocTime, eMallocTime;
    TIMESTAMP( sMallocTime );
#endif
    struct Params inputParam;
    /// Get the field information.
    nfiields = mxGetNumberOfFields(prhs[0]);
    //TODO: add an assertion to check the value of nfields
    maxElem = 2500; /// empirical value and is subject to change.
    
    sizeTol = sizeLim = 0;
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
    
    /// Assignment of the struct inputParam.   
    inputParam.Iph = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,0)));
    inputParam.Io = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,1)));
    inputParam.Rs = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,2)));
    inputParam.nVth = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,3)));
    inputParam.Rsh = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,4)));
    inputParam.di2Mutau = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,5)));
    inputParam.Vbi = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,6)));
    inputParam.Brev = *((double*)mxGetData(mxGetFieldByNumber(prhs[0],0,7)));
    
    /// Output
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
    
    mexODMpwlapprox(inputParam, tol, lim, sizeTol, sizeLim, outMatrix1, outMatrix2, *(unsigned int*)mxGetData(plhs[2]));
#ifdef TIMER
    TIMESTAMP( eCallTime );
    std::cout<<"Call overhead...";
    std::cout<<"Elapsed time is "<< eCallTime - sCallTime <<" seconds\n"; 
#endif

#ifdef TIMER
    TIMESTAMP( eWholeTime );
    std::cout<<"Whole overhead...";
    std::cout<<"Elapsed time is "<< eWholeTime - sWholeTime <<" seconds\n";
#endif
}
