#!/bin/bash
echo "cal gofr"
./gofr.out ../BCC.xyz ./gofrBCC.txt 1
echo "BCC Cry OK"
./gofr.out ../HCP.xyz ./gofrHCP.txt 1
echo "HCP Cry OK"
./gofr.out ../FCC.xyz ./gofrFCC.txt 1
echo "FCC Cry OK"
./gofr.out ../P0.01T0.006.xyz ./gofrP0.01T0.006.txt 1000
./gofr.out ../P1T0.005.xyz ./gofrP1T0.005.txt 1000
echo "GAUSS Liq OK"
./gofr.out ../P1T0.0018.xyz ./gofrP1T0.0018.txt 1000
./gofr.out ../P0.05T0.0052.xyz ./gofrP0.05T0.0052.txt 1000
echo "GAUSS BCC OK"
./gofr.out ../P0.01T0.00262.xyz ./gofrP0.01T0.00262.txt 1000
echo "GAUSS FCC OK"
./gofr.out ../P1T0.1.xyz ./gofrP1T0.1.txt 1000
echo "LJ HCP OK"
exit 0
