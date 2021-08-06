#!/bin/bash
echo "cal conf"
for i in $(seq 1 9)
do
	./gaussian.out Rbcc.xyz bcc 0.0"$i" bcc0"$i".txt
	./gaussian.out Rfcc.xyz fcc 0.0"$i" fcc0"$i".txt
	./gaussian.out Rhcp.xyz hcp 0.0"$i" hcp0"$i".txt
	echo "R"$i" cal ok"
done
for i in $(seq 10 20)
do
	./gaussian.out Rbcc.xyz bcc 0."$i" bcc"$i".txt
	./gaussian.out Rfcc.xyz fcc 0."$i" fcc"$i".txt
	./gaussian.out Rhcp.xyz hcp 0."$i" hcp"$i".txt
	echo "R"$i" cal ok"
done
exit 0




