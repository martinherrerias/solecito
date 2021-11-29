#include "mex.hpp"
#include "mexAdapter.hpp"
#include "../pvProj/ElectricalCalculation/MonoDescPWL.h"
#include "../pvProj/testUtils.h"

//using namespace Eigen;
class MonoDescPWL;

using matlab::mex::ArgumentList;
using namespace matlab::data;


class MexFunction : public matlab::mex::Function {
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        int i, j;
        Eigen::MatrixXd PS, P2;
        
        GetSimOptFile* filePtr = GetSimOptFile::Instance();
        filePtr->openOptionFile("../pvProj/Resources/SimOption.xml"); //TODO: the path should be changed
        filePtr->readOptionFile();
        
        // Move object to variable
        matlab::data::Array blockArray = std::move(inputs[0]);
        matlab::data::Array diodeArray = std::move(inputs[1]);
        
        ArrayDimensions inputSize = blockArray.getDimensions();
        int rows = inputSize[0];
        int cols = inputSize[1];
        
        std::vector<MonoDescPWL*> set;
        std::vector<std::vector<MonoDescPWL*>> result(inputSize[1], std::vector<MonoDescPWL*>(inputSize[0]));
        MonoDescPWL pwl1, pwl2;
        std::vector<double> bound;
  
        double sModuleTime, eModuleTime;
        TIMESTAMP( sModuleTime );
        
        double sInterTime, eInterTime, getPropertySum = 0.0;
        double sAssignTime, eAssignTime, assignSum = 0.0;
        
        for (j = 0; j < rows*cols; j++)
       {
            TIMESTAMP( sInterTime ); 
            const matlab::data::Array voltageBlock = matlabPtr->getProperty(blockArray, j, u"x");
            const matlab::data::TypedArray<double> currentBlock = matlabPtr->getProperty(blockArray, j, u"y");

            ArrayDimensions blockSize = voltageBlock.getDimensions();

            const matlab::data::TypedArray<double> voltageDiode = matlabPtr->getProperty(diodeArray, j, u"x");
            const matlab::data::TypedArray<double> currentDiode = matlabPtr->getProperty(diodeArray, j, u"y");

            TIMESTAMP( eInterTime );
            getPropertySum += eInterTime - sInterTime;

            matlab::data::ArrayDimensions diodeSize = voltageDiode.getDimensions();

            PS.resize(blockSize[0],2);
            P2.resize(diodeSize[0],2);

            TIMESTAMP( sAssignTime);
            for (i = 0; i < blockSize[0]; i++)
            {
                PS(i, 0) = voltageBlock[i];
                PS(i, 1) = currentBlock[i];
            }   
            for(i = 0; i < diodeSize[0]; i++)
            {
                P2(i, 0) = voltageDiode[i];
                P2(i, 1) = currentDiode[i];
            }

            TIMESTAMP(eAssignTime);
            assignSum += eAssignTime - sAssignTime;

            pwl1 = MonoDescPWL(PS);
            pwl2 = MonoDescPWL(P2);

            set.push_back(&pwl1);
            set.push_back(&pwl2);

             /* if (j == 0)
            {
                std::cout<<"pwl1:\n";
                std::cout<<pwl1.pointSet<<"\n";
                std::cout<<"pwl2:\n";
                std::cout<<pwl2.pointSet<<"\n";
            }*/
            result[j/rows][j%rows] = new MonoDescPWL(set, Dim::vert, {EPSILON}, {-Inf, Inf, 18.1462,  -9.0731});

            /*  if (j == 0)
            {
                std::cout<<"result is \n";
                std::cout<<result[0][0]->pointSet<<"\n";
            }*/
            result[j/rows][j%rows]->getMpp(bound);
            set.clear();
       }
       TIMESTAMP( eModuleTime );
       std::cout<<"Adding bypass diodes and getting isolated MPPs...";
       std::cout<<"the summe is "<<getPropertySum<<"\n";
       std::cout<<"the sum is "<<assignSum<<"\n";
       std::cout<<"The elapsed time is "<< eModuleTime - sModuleTime << " seconds\n";
        
        PS.resize(0,0);
        P2.resize(0,0);
    }
};
