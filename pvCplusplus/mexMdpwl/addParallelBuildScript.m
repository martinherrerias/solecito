% addparallel without TIMER defined
mex -c ./mexaddparallel.cpp CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../EIGENDIR -I../LIBIGLDIR
mex CFLAGS="\$CFLAGS -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I. mexaddparallel.o MonoDescPWL.o xmlutils.o

%%
% addparallel with TIME defined
%{
mex -c ./mexaddparallel.cpp CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I../EIGENDIR 
mex CFLAGS="\$CFLAGS -DTIMER -O3 -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"  -I. mexaddparallel.o MonoDescPWL.o xmlutils.o
%}
