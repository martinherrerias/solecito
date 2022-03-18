#include "../pvProj/ElectricalCalculation/MonoDescPWL.h"
#include "../pvProj/testUtils.h"
#include "mex.h"
#include <omp.h>
#include <sched.h>
#include <string.h>
#include <limits.h> /* PATH_MAX */


using namespace Eigen;
class MonoDescPWL;

void mexaddparallel(double ***array, mwSize* inLength, mwSize size, double* tol, double* lim, mwSize tolSize, mwSize limSize, double* outarray1, double* outarray2, unsigned int& length)
{
#ifdef TIMER
    double sXMLTime, eXMLTime;
    TIMESTAMP( sXMLTime );
#endif
    /// Get file pointer to the SimOption.xml file.
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

    double sModuleTime, eModuleTime;
    double sIntTime, eIntTime;
    TIMESTAMP( sModuleTime );
#endif
    MatrixXd PS;
    int i, j;
    
    std::vector<MonoDescPWL> pwlset((size));
    std::vector<MonoDescPWL*> set((size));
    std::vector<double> vecTol;
    std::vector<double> vecLim;
    MonoDescPWL output;
       
    if (tolSize != 0)
        vecTol.assign(tol, tol+tolSize);
    if (limSize != 0)
        vecLim.assign(lim, lim+limSize);
  
#pragma omp parallel private(j, PS) if (size > 32)
    {
#pragma omp for
    for (j = 0; j < size; j++)
    {
        PS.resize(inLength[j],2);     
        for(i = 0; i < inLength[j]; i++)
        {
            PS(i, 0) = array[j][0][i];
            PS(i, 1) = array[j][1][i];
        }
        /// Create a mdpwl obj based on each array.
        pwlset[j] = MonoDescPWL(PS, {0.0});
        set[j] = &pwlset[j];
    }
    }
#ifdef TIMER
    TIMESTAMP( eModuleTime );
#endif
    /// addparallel
    output = MonoDescPWL(set, Dim::vert, vecTol, vecLim);
 #ifdef TIMER
    std::cout<<"Pure C++ addparallel overhead...";
    std::cout<<"Elapsed time is "<< eModuleTime - sModuleTime <<" seconds\n";

    double sOutputTime, eOutputTime;
    TIMESTAMP( sOutputTime );
#endif
    length = output.pointSet.rows();
    /// Extract the array from the created mdpwl obj.
    for (j = 0; j < length; j++)
    {
        outarray1[j] = output.pointSet(j,0);
        outarray2[j] = output.pointSet(j,1);
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
    double *outMatrix1;
    double *outMatrix2;
    double ***inMatrix, *tol, *lim;
    mwSize *inLength;

    mwSize i, j, nrows0, ncols0, nrows1, ncols1, size, maxElem, sizeTol, sizeLim;
#ifdef TIMER
    double sMallocTime, eMallocTime;
    TIMESTAMP( sMallocTime );
#endif
    size = mxGetScalar(prhs[0]); /// The number of the passed arrays
    inMatrix = (double***)mxMalloc( size* sizeof(double**) );
    inLength = (mwSize*)mxMalloc( size*sizeof(mwSize) );
    
    /// Fill tol and lim if there are in the prhs.
    sizeTol = sizeLim = 0;
    if (nrhs >= size+1)
    {
        tol = mxGetPr(prhs[size+1]);
        sizeTol = mxGetN(prhs[size+1]);
    }
    if (nrhs >= size+2)
    {
        lim = mxGetPr(prhs[size+2]);
        sizeLim = mxGetN(prhs[size+2]);
    }
    
    maxElem = 0;
    
    for(i = 1; i <= size; i++)
    {
        nrows0 = mxGetM(prhs[i]);
        ncols0  = mxGetN(prhs[i]);
        inMatrix[i - 1] = (double**) mxMalloc( ncols0* sizeof(double*) );
        inMatrix[i - 1][0] = mxGetPr(prhs[i]);
        for (j = 1; j< ncols0; j++)
            inMatrix[i - 1][j] = inMatrix[i - 1][j-1] + nrows0;
        inLength[i-1] = nrows0;
        maxElem += nrows0;
    }
    
//    tol = mxGetPr(prhs[size+1]);
//    lim = mxGetPr(prhs[size+2]);
    
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
    mexaddparallel(inMatrix, inLength,  size, tol, lim, sizeTol, sizeLim, outMatrix1, outMatrix2, *(unsigned int*)mxGetData(plhs[2]));
#ifdef TIMER
    TIMESTAMP( eCallTime );
    std::cout<<"call overhead...";
    std::cout<<"Elapsed time is "<< eCallTime - sCallTime <<" seconds\n"; 
#endif
    for (i = 0; i < size; i++)
        mxFree (inMatrix[i]);
    mxFree (inMatrix);
#ifdef TIMER
    TIMESTAMP( eWholeTime );
    std::cout<<"whole overhead...";
    std::cout<<"Elapsed time is "<< eWholeTime - sWholeTime <<" seconds\n";
#endif
}
