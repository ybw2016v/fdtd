fdtdcal.so:fdtdcal.c
	gcc fdtdcal.c -lm -fopenmp -fPIC -shared -o fdtdcal.so