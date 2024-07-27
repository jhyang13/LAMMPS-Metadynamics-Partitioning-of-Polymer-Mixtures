#!/bin/bash 

NP=4

for NA in 100
do
	for PHIA in 1
	do
		for NB in 20
		do	
			for PHIB in 1
			do
				for R in 6
				do
					for D in 12
					do
						for EPS in 0.40 0.60 0.80
						do 	
							for DATE in 20230509
							do
					
									DIR="${DATE}_NA${NA}_PHIA${PHIA}_NB${NB}_PHIB${PHIB}_R${R}_EPS${EPS}_hill0.01"
									mkdir $DIR
									cp templates/zcoord.colvars templates/def.chain templates/poly.data templates/in.mono templates/submit.mino $DIR
									cd $DIR
									
									RAND1=$RANDOM
 		                                                        RAND2=$RANDOM

									sed '{s/{EPS}/'"$EPS"'/g}' < in.mono > tmp1
									sed '{s/{D}/'"$D"'/g}' < tmp1 > tmp2
									sed '{s/{R}/'"$R"'/g}' < tmp2 > tmp3
									sed '{s/{RAND1}/'"$RAND1"'/g}' < tmp3 > tmp4
									sed '{s/{RAND2}/'"$RAND2"'/g}' < tmp4 > tmp5

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
			done
		done
	done
done
