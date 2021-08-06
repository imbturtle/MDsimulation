#!/bin/bash
echo "cal nei"
./neiborlist.out ../../P0.01T0.006.xyz GAUSSLiq NB0P0.01T0.006.txt CN0P0.01T0.006.txt 200 16 1.5







'''
./neiborlist.out ../../P0.01T0.006.xyz GAUSS NB0P0.01T0.006.txt CN0P0.01T0.006.txt 500 16 1.5
./neiborlist.out ../../P1T0.005.xyz GAUSS NB0P1T0.005.txt CN0P1T0.005.txt 500 16 1.5
echo "GAUSS Liq OK"
./neiborlist.out ../../P1T0.0018.xyz GAUSS NB1P1T0.0018.txt CN1P1T0.0018.txt 500 16 1.5
./neiborlist.out ../../P0.05T0.0052.xyz GAUSS NB1P0.05T0.0052.txt CN1P0.05T0.0052.txt 500 16 1.5
echo "GAUSS BCC OK"
./neiborlist.out ../../P0.01T0.00262.xyz GAUSS NB2P0.01T0.00262.txt CN2P0.01T0.00262.txt 500 16 1.5
echo "GAUSS FCC OK"
./neiborlist.out ../../P1T0.1.xyz LJ NB3P1T0.1.txt CN3P1T0.1.txt 500 16 1.5
echo "LJ HCP OK"


exit 0

#!/bin/bash
echo "cal gofr"
./gofr.out ../BCC.xyz LJ ./gofrBCC.txt 1
echo "BCC Cry OK"
./gofr.out ../HCP.xyz LJ ./gofrHCP.txt 1
echo "HCP Cry OK"
./gofr.out ../FCC.xyz LJ ./gofrFCC.txt 1
echo "FCC Cry OK"
./gofr.out ../P0.01T0.006.xyz GAUSS ./gofrP0.01T0.006.txt 1000
./gofr.out ../P1T0.005.xyz GAUSS ./gofrP1T0.005.txt 1000
echo "GAUSS Liq OK"
./gofr.out ../P1T0.0018.xyz GAUSS ./gofrP1T0.0018.txt 1000
./gofr.out ../P0.05T0.0052.xyz GAUSS ./gofrP0.05T0.0052.txt 1000
echo "GAUSS BCC OK"
./gofr.out ../P0.01T0.00262.xyz GAUSS ./gofrP0.01T0.00262.txt 1000
echo "GAUSS FCC OK"
./gofr.out ../P1T0.1.xyz LJ ./gofrP1T0.1.txt 1000
echo "LJ HCP OK"
exit 0
'''
