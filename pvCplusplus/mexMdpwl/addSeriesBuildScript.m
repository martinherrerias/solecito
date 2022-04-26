% addseries without TIMER defined
mex -c ./mexaddseries.cpp CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../eigen -I../libigl
mex CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS  -fopenmp"  -I. mexaddseries.o MonoDescPWL.o xmlutils.o

%% 
% addseries with TIMER defined
%{
mex -c ./mexaddseries.cpp CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../eigen 
mex CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS  -fopenmp"  -I. mexaddseries.o MonoDescPWL.o xmlutils.o
%}
