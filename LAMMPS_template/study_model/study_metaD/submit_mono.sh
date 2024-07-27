#!/bin/bash 

NP=4

for N in 120
do	
	for R in 3
	do
		for D in 6
		do
			for EPS in 0.60
			do 
				for DATE in 20230105
				do
					for REP in 01 02 03 04 05
					do
									
						DIR="${DATE}_N${N}_R${R}_EPS${EPS}_REP${REP}"
						mkdir $DIR
						cp templates/def.chain templates/poly.data templates/in.mono templates/submit.mino templates/zcoord.colvars $DIR
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
						
						sed '{s/{N}/'"$N"'/g}' < zcoord.colvars > tmp1
						
						mv tmp1 zcoord.colvars
						rm tmp*
 					
						qsub submit.mino
						cd ..
					done
	 			done			
			done
		done
	done
done
