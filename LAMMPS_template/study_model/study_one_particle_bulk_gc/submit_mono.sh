#!/bin/bash 

NP=1

for NA in 1
do
	for u in -3.11 -3.12 -3.13 -3.14 
	do 
		for EPS in 0.00
		do 
			for DATE in 20230129
			do
									
				DIR="${DATE}_NA${NA}_EPS${EPS}_u${u}"
				mkdir $DIR
				cp templates/def.chain templates/minsert.txt templates/poly.data templates/in.mono templates/submit.mino $DIR
				cd $DIR
									
				RAND1=$RANDOM
 		        	RAND2=$RANDOM
				RAND3=$RANDOM
			
				sed '{s/{EPS}/'"$EPS"'/g}' < in.mono > tmp1
				sed '{s/{RAND1}/'"$RAND1"'/g}' < tmp1 > tmp2
				sed '{s/{RAND2}/'"$RAND2"'/g}' < tmp2 > tmp3
				sed '{s/{RAND3}/'"$RAND3"'/g}' < tmp3 > tmp4
				sed '{s/{u}/'"$u"'/g}' < tmp4 > tmp5

				mv tmp5 in.run
				rm tmp* in.mono

				sed '{s/{NP}/'"$NP"'/g}' < submit.mino > tmp1
				sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
						
				mv tmp2 submit.mino
				rm tmp*	
											
				qsub submit.mino
				cd ..
			done
		done
	done
done
