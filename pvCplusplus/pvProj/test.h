#ifndef TEST_H_INCLUDED
#define TEST_H_INCLUDED
#define PRECISION 0.0001

#define EXPECTMATRIX_EQ( arg1_, arg2_ ) (((arg1_)-(arg2_)).norm() < PRECISION)
#define EXPECTNUM_EQ( arg1_, arg2_ ) (((arg1_) - (arg2_)) < PRECISION)

#endif
