% addseries without TIMER defined
mex -c ./mexODMpwlapprox.cpp CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../eigen -I../libigl
mex -c ./mexVecODMpwlapproxAddseries.cpp CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../eigen -I../libigl
mex CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS  -fopenmp"  -I.  mexODMpwlapprox.o OneDiodeModel.o xmlutils.o LambertW.o
mex CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS  -fopenmp"  -I.  mexVecODMpwlapproxAddseries.o OneDiodeModel.o MonoDescPWL.o xmlutils.o LambertW.o

%% 
% addseries with TIMER defined
%{
mex -c ./mexaddseries.cpp CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../eigen 
mex CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS  -fopenmp"  -I. mexaddseries.o MonoDescPWL.o xmlutils.o
%}
