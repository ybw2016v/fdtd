fdtdcal.so:fdtdcal.c fdtdcal.h calmur2.c
	gcc calmur2.c fdtdcal.c -lm -fopenmp -fPIC -shared -O3 -o fdtdcal.so
upload:
	scp fdtdcal.c fdtdcal.h txx.py fdtdcallib.py plot.py calmur2.c Makefile guest@jn.dogcraft.top:/home/guest/fdog/
download:
	scp guest@jn.dogcraft.top:/home/guest/fdog/*.jpg ./
