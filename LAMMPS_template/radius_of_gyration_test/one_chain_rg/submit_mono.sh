#!/bin/bash 

NP=1

for N in 2
do	 
	for DATE in 20230406
	do
									
		DIR="${DATE}_N${N}_rg"
		mkdir $DIR
		cp templates/def.chain templates/poly.data templates/in.mono templates/submit.mino $DIR
		cd $DIR
									
		RAND1=$RANDOM
 		RAND2=$RANDOM

		sed '{s/{RAND1}/'"$RAND1"'/g}' < in.mono > tmp1
		sed '{s/{RAND2}/'"$RAND2"'/g}' < tmp1 > tmp2

		mv tmp2 in.run
		rm tmp* in.mono

		sed '{s/{NP}/'"$NP"'/g}' < submit.mino > tmp1
		sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
		mv tmp2 submit.mino
		rm tmp*	
						
		qsub submit.mino
		cd ..

	done
done
