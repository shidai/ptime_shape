#!/bin/sh

gcc -Wall -lm -lfftw3 -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptime_shape.c readfits.c shape.c taylor92.c -o ptime_shape.out
#gcc -Wall -lm -lfftw3 -lgsl -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio ptime_shape.c readfits.c shape.c taylor92.c -o ptime_shape.out
