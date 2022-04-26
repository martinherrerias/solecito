%Compile the common object files without TIMER defined
mex -c ../pvProj/Utils/xmlutils.cpp CFLAGS="\$CFLAGS -O3 -fopenmp"  LDFLAGS="\$LDFLAGS -fopenmp" -I../eigen -I../libigl -I../pvProj/Utils
mex -c ../pvProj/ElectricalCalculation/MonoDescPWL.cpp  CFLAGS="\$CFLAGS -O3 -fopenmp"  LDFLAGS="\$LDFLAGS -fopenmp" -I../eigen -I../libigl -I../pvProj/Utils
mex -c ../pvProj/Utils/LambertW.cpp CFLAGS="\$CFLAGS -O3 -fopenmp"  LDFLAGS="\$LDFLAGS -fopenmp" -I../eigen -I../libigl -I../pvProj/Utils
mex -c ../pvProj/ModuleDefinition/OneDiodeModel.cpp  CFLAGS="\$CFLAGS -O3 -fopenmp"  LDFLAGS="\$LDFLAGS -fopenmp" -I../eigen -I../libigl -I../pvProj/Utils
