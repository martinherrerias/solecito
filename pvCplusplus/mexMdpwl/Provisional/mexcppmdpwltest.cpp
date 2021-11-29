#include "../pvProj/ElectricalCalculation/MonoDescPWL.h"
#include "../pvProj/testUtils.h"
#include "mex.h"
#include <omp.h>
#include <sched.h>
#include <string.h>
#define MaxElem 2500

using namespace Eigen;
class MonoDescPWL;

/*static char *cpuset_to_cstr(cpu_set_t *mask, char *str)
{
  char *ptr = str;
  int i, j, entry_made = 0;
  for (i = 0; i < CPU_SETSIZE; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < CPU_SETSIZE; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d,", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d,", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d,", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  *ptr = 0;
  return(str);
}*/

void mexcppmdpwltest(double **array0, double**array1, double* outarray1, double* outarray2, double e, double* m0, double* m1)
{// inside the mexmdpwl, we perform the addseries or addparallel
    GetSimOptFile* filePtr = GetSimOptFile::Instance();
    filePtr->openOptionFile("../pvProj/Resources/SimOption.xml");
    filePtr->readOptionFile();

    int lengthupm0, lengthupm1;

    MatrixXd PS, P2;
    int i, j;
    int rows, cols;
    rows = 20;
    cols = 6;
    
    std::vector<MonoDescPWL*> set;
    std::vector<std::vector<MonoDescPWL*>> result(cols, std::vector<MonoDescPWL*>(rows));
    std::vector<MonoDescPWL*> medium(cols);
    MonoDescPWL pwl1, pwl2, output;
    std::vector<double> bound;

    int end = rows*cols; //rows*cols
    std::vector<int> lengthm0(end), lengthm1(end);
    
    lengthm0[0] = (int)m0[0];
    lengthm1[0] = (int)m1[0];
    for (j = 1; j < end; j++)
    {
        lengthm0[j] = lengthm0[j - 1] + (int)m0[j];
        lengthm1[j] = lengthm1[j - 1] + (int)m1[j];
    }
    
    double sModuleTime, eModuleTime;
    double sIntTime, eIntTime;
 
    TIMESTAMP( sModuleTime );
#pragma omp parallel private(j, lengthupm0, lengthupm1, i, PS, P2, pwl1, pwl2, set)
    {
#pragma omp for
    for (j = 0; j < end; j++)
    {   
        lengthupm0 = (int)m0[j];
        lengthupm1 = (int)m1[j];
        PS.resize(lengthupm0,2);
        P2.resize(lengthupm1,2);
        
        int k = 0;
        if (j == 0)
            i = 0;
        else i = lengthm0[j-1];
        for (; i < lengthm0[j]; i++)
        {
            PS(k, 0) = array0[0][i];
            PS(k, 1) = array0[1][i];
            k++;
        }   
                
        k = 0;
        if (j == 0)
            i = 0;
        else i = lengthm1[j-1];
        for (; i < lengthm1[j]; i++)
        {
            P2(k, 0) = array1[0][i];
            P2(k, 1) = array1[1][i];
            k++;
        }
    
        pwl1 = MonoDescPWL(PS);
        pwl2 = MonoDescPWL(P2);

        set.push_back(&pwl1);
	    set.push_back(&pwl2);
        
  /*       if (j == 0)
        {
            std::cout<<"pwl1:\n";
            std::cout<<pwl1.pointSet<<"\n";
            std::cout<<"pwl2:\n";
            std::cout<<pwl2.pointSet<<"\n";
        }*/
        
        result[j/rows][j%rows] = new MonoDescPWL(set, Dim::vert, {EPSILON0}, {-Inf, Inf, 18.1462,  -9.0731});
        
/*        if (j == 0)
        {
            std::cout<<"result is \n";
            std::cout<<result[0][0]->pointSet<<"\n";
        }*/
        result[j/rows][j%rows]->getMpp(bound);
        
        set.clear();
    }
    }
    TIMESTAMP( eModuleTime );
    std::cout<<"Adding bypass diodes and getting isolated MPPs...";
    std::cout<<"Elapsed time is "<< eModuleTime - sModuleTime << " seconds\n";
    
//   std::cout<<"new 0:\n"<<result[0][0]->pointSet<<"\n";
//   std::cout<<"new 1:\n"<<result[0][1]->pointSet<<"\n";
  
    double sStringTime, eStringTime;
    TIMESTAMP( sStringTime );
    
/*    int thread;
    int numthreads;
    cpu_set_t coremask;
    char clbuf[7 * CPU_SETSIZE], hnbuf[64];
    
    memset(clbuf, 0, sizeof(clbuf));
    memset(hnbuf, 0, sizeof(hnbuf));
    (void)gethostname(hnbuf, sizeof(hnbuf));*/
    
#pragma omp parallel private(j)
    {
    /*    thread = omp_get_thread_num();
        numthreads = omp_get_num_threads();
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, clbuf);
        #pragma omp barrier
        
        #pragma omp critical
        std::cout<<"Hello from thread "<<thread<<". (core affinity = "<<clbuf<<" )"<<", the num of threads is "<<numthreads<<"\n";
    */
#pragma omp for
    for (j = 0; j < cols; j++)
    {
        medium[j] = new MonoDescPWL(result[j], Dim::hori);
    }    
    }
//    std::cout<<"hori0:\n"<<medium[0]->pointSet<<"\n";
    
    TIMESTAMP( eStringTime );

    std::cout<<"Adding blocks in series...";
    std::cout<<"Elapsed time is "<< eStringTime - sStringTime << " seconds\n";
    
    double sArrayTime, eArrayTime;
    TIMESTAMP( sArrayTime );
    output = MonoDescPWL(medium, Dim::vert, {EPSILON0}, {-Inf, Inf, 181.4622, 0});
    output.getMpp(bound);
    TIMESTAMP( eArrayTime );
    
    std::cout<<"Adding strings in parallel and getting array MPP...";
    std::cout<<"Elapsed time is "<< eArrayTime - sArrayTime << " seconds\n";
    
    double sOutputTime, eOutputTime;
    TIMESTAMP( sOutputTime );

    for (j = 0; j < output.pointSet.rows(); j++)
    {
        outarray1[j] = output.pointSet(j,0);
        outarray2[j] = output.pointSet(j,1);
    }

    TIMESTAMP( eOutputTime );
    std::cout<<"Output overhead...";
    std::cout<<"Elapsed time is "<< eOutputTime - sOutputTime <<" seconds\n";     
}
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *outMatrix1;
    double *outMatrix2;
    double **inMatrix0, **inMatrix1;
    double *input;
    double e = 1.0;

    mwSize i, j, nrows0, ncols0, nrows1, ncols1, index;
    double *num0, *num1;

    nrows0 = mxGetM(prhs[0]);
    ncols0  = mxGetN(prhs[0]);
    
    nrows1 = mxGetM(prhs[1]);
    ncols1  = mxGetN(prhs[1]);
    
    index    = mxGetN(prhs[2]);

    inMatrix0 = (double**) mxMalloc( ncols0* sizeof(double*) );
    inMatrix0[0] = mxGetPr(prhs[0]);
    for (j = 1; j< ncols0; j++)
        inMatrix0[j] = inMatrix0[j-1] + nrows0;
    
    inMatrix1 = (double**) mxMalloc( ncols1* sizeof(double*) );
    inMatrix1[0] = mxGetPr(prhs[1]);
    for (j = 1; j< ncols1; j++)
        inMatrix1[j] = inMatrix1[j-1] + nrows1;
    
    num0 = (double*)mxGetPr(prhs[2]);
    num1 = (double*)mxGetPr(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(1, MaxElem, mxREAL);
    outMatrix1 = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, MaxElem, mxREAL);
    outMatrix2 = mxGetPr(plhs[1]);
    mexcppmdpwltest(inMatrix0, inMatrix1,  outMatrix1, outMatrix2, e, num0, num1);
    mxFree (inMatrix0);
    mxFree (inMatrix1);
}
