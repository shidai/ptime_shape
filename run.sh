#!/bin/sh

gcc -lm -lfftw3 -lgsl -L/usr/local/lib/cfitsio -I/usr/include/cfitsio/ -lcfitsio main.c readfits.c shape.c taylor92.c -o ptime_shape.out
