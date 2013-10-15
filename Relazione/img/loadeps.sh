#!/bin/bash

#Creo le liste di cartelle  che voglio compilare
	find ./* -type d | sed 's/.\///' > cartelle.txt
	#Determino il numero di file
	NFILE=`find ./* -type d | wc -l`
	for((i=1;i<=$NFILE;i++));do 
	cd `sed -n "$i"p cartelle.txt`;
	echo "`sed -n "$i"p ../cartelle.txt`";
	cd ../; done
