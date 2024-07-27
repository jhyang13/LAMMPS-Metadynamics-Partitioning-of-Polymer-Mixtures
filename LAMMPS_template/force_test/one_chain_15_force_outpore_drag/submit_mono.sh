#!/bin/bash 

NP=1

for N in 15
do	
	for R in 2
	do
		for D in 4
		do
			for EPS in 0.20
			do 
				for DATE in 20230331
				do
					#for h in -20 -15 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0
					#do
					#for h in 1 2 3 4 5 6 7 8 9 10 15 20
                                        #do
					#for h in -10
                                        #do
									
						DIR="${DATE}_N${N}_R${R}_EPS${EPS}_lad10"
						mkdir $DIR
						cp templates/poly.data templates/in.mono templates/submit.mino $DIR
						cd $DIR
									
						RAND1=$RANDOM
 		                                RAND2=$RANDOM

						sed '{s/{EPS}/'"$EPS"'/g}' < in.mono > tmp1
						sed '{s/{D}/'"$D"'/g}' < tmp1 > tmp2
						sed '{s/{R}/'"$R"'/g}' < tmp2 > tmp3
						sed '{s/{RAND1}/'"$RAND1"'/g}' < tmp3 > tmp4
						sed '{s/{RAND2}/'"$RAND2"'/g}' < tmp4 > tmp5
						sed '{s/{h}/'"$h"'/g}' < tmp5 > tmp6

						mv tmp6 in.run
						rm tmp* in.mono

						sed '{s/{NP}/'"$NP"'/g}' < submit.mino > tmp1
						sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
						mv tmp2 submit.mino
						rm tmp*	
						
						#qsub submit.mino
						cd ..

	 			done			
			done
		done
	done
done
