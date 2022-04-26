%Compiler flag for mex object
mex  -v CFLAGS="\$CFLAGS -funroll-loops -O3"  ./mexcppmdpwltest_obj.cpp ../pvProj/Utils/xmlutils.cpp ../pvProj/ElectricalCalculation/MonoDescPWL.cpp -I.  -I../pvProj/Utils -I../eigen 

