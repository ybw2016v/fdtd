fdtdcal.so:fdtdcal.c fdtdcal.h
	gcc fdtdcal.c -lm -fopenmp -fPIC -shared -o fdtdcal.so