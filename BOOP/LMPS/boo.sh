#!/bin/bash
echo "cal boo"					
./boo.out ../BCC.xyz LJ 1.135 ./BOOPBCC.txt 1
echo "BCC Cry OK"
./boo.out ../HCP.xyz LJ 1.165 ./BOOPHCP.txt 1
echo "HCP Cry OK"
./boo.out ../FCC.xyz LJ 1.165 ./BOOPFCC.txt 1
echo "FCC Cry OK"
./boo.out ../P0.01T0.006.xyz GAUSS 3.265 ./BOOPP0.01T0.006.txt 1
./boo.out ../P1T0.005.xyz GAUSS 1.845 ./BOOPP1T0.005.txt 1
echo "GAUSS Liq OK"
./boo.out ../P1T0.0018.xyz GAUSS 1.815 ./BOOPP1T0.0018.txt 1
./boo.out ../P0.05T0.0052.xyz GAUSS 2.755 ./BOOPP0.05T0.0052.txt 1
echo "GAUSS BCC OK"
./boo.out ../P0.01T0.00262.xyz GAUSS 2.895 ./BOOPP0.01T0.00262.txt 1
echo "GAUSS FCC OK"
./boo.out ../P1T0.1.xyz LJ 1.245 ./BOOPP1T0.1.txt 1
echo "LJ HCP OK"
exit 0

