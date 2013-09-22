#!/bin/bash

echo "Cosa vuoi fare?(CREA/RIMUOVI/CONTA)"
read ORDINE

CREA="CREA"
RIMUOVI="RIMUOVI"
CONTA="CONTA"

if [ "$ORDINE" == "$CREA" ]; then
	#Creo le liste di file che voglio convertire
	ls *.png > files_png.txt
	sed 's/.png/.eps/' files_png.txt > files_eps.txt
	sed 's/.png//' files_png.txt > files.txt
	#Determino il numero di file
	NFILE=`find *.png -type f | wc -l`
	for((i=1;i<=$NFILE;i++));do convert `sed -n "$i"p files_png.txt` `sed -n "$i"p files_eps.txt`; done
	git add .
else if [ "$ORDINE" == "$RIMUOVI" ]; then
	git rm -f files*
	git rm -f *.eps
else if [ "$ORDINE" == "$CONTA" ]; then
	echo "Numero file png"
	find *.png -type f | wc -l
	echo "Numero file eps"
	find *.eps -type f | wc -l
	echo "Numero file txt"
	find *.txt -type f | wc -l
fi
fi
fi
