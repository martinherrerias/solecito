%Compiler flag supporting openmp
mex -c ./mexcppmdpwltest.cpp CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../EIGENDIR 
mex  CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I. mexcppmdpwltest.o MonoDescPWL.o xmlutils.o

%Compiler flag ignoring openmp
%mex -c ../pvProj/Utils/xmlutils.cpp CFLAGS="\$CFLAGS -DNDEBUG -funroll-loops  -ffast-math -O3"  -I../EIGENDIR -I../pvProj/Utils
%mex -c ../pvProj/ElectricalCalculation/MonoDescPWL.cpp  CFLAGS="\$CFLAGS -DNDEBUG -funroll-loops -ffast-math -O3" -I../EIGENDIR -I../pvProj/Utils
%mex  CFLAGS="\$CFLAGS -DNDEBUG -funroll-loops -ffast-math -O3"  -I../EIGENDIR -c ./mexcppmdpwltest.cpp
%mex  CFLAGS="\$CFLAGS -DNDEBUG -funroll-loops -ffast-math -O3"  -I. mexcppmdpwltest.o MonoDescPWL.o xmlutils.o

%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -I../pvProj/Utils -c ../pvProj/Utils/xmlutils.cpp
%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -I../pvProj/Utils -c ../pvProj/ElectricalCalculation/MonoDescPWL.cpp

%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -c ./mexaddparallel.cpp
%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -c ./mexaddseries.cpp
%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -I../pvProj/Utils -c ./mexcppmdpwltest.cpp

%mex COMPFLAGS='$COMPFLAGS -DNDEBUG -funroll-loops -march=native -ffast-math -O3' -I../EIGENDIR -c ./mexaddparallelobj.cpp

%mex COMPFLAGS='$COMPFLAGS -DNDEBUG  -funroll-loops -march=native -ffast-math -O3' -I. mexaddparallel.o  MonoDescPWL.o xmlutils.o 
%mex COMPFLAGS='$COMPFLAGS -DNDEBUG  -funroll-loops -march=native -ffast-math -O3' -I. mexaddseries.o MonoDescPWL.o xmlutils.o
