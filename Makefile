fdtdcal.so:fdtdcal.c fdtdcal.h
	gcc fdtdcal.c -lm -fopenmp -fPIC -shared -O3 -o fdtdcal.so
upload:
	scp fdtdcal.c fdtdcal.h txx.py fdtdcallib.py plot.py guest@jn.dogcraft.top:/home/guest/fdog/
download:
	scp guest@jn.dogcraft.top:/home/guest/fdog/*.jpg ./
