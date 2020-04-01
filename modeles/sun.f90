#!/bin/sh
#-Xlist

LIB=/users-data/chauvin/SIMBOLX/F90/lib

rm $1
f90 -V -C -O3 -M$LIB  $1.f90 -o $1 /opt/NAGnl/flsol19da/libnag.a
rm $1.o

 
