#!/bin/bash
echo "cal bad"
./bad.out ../BCC.xyz LJ 1.295 ./BadBCC.txt ./CNBCC.txt 1
echo "BCC Cry OK"
./bad.out ../HCP.xyz LJ 1.165 ./BadHCP.txt ./CNHCP.txt 1
echo "HCP Cry OK"
./bad.out ../FCC.xyz LJ 1.165 ./BadFCC.txt ./CNFCC.txt 1
echo "FCC Cry OK"
./bad.out ../P1T0.1.xyz LJ 1.275 ./BadP1T0.1.txt ./CNP1T0.1.txt 1000
echo "LJ HCP OK"
./bad.out ../P1T0.0018.xyz LJ 1.815 ./BadP1T0.0018.txt ./badP1T0.0018.txt 1000
./bad.out ../P0.05T0.0052.xyz LJ 2.755 ./BadP0.05T0.0052.txt ./CNP0.05T0.0052.txt 1000
echo "GAUSS BCC OK"
./bad.out ../P0.01T0.00262.xyz LJ 2.915 ./BadP0.01T0.00262.txt ./CNP0.01T0.00262.txt 1000
echo "GAUSS FCC OK"
exit 0

