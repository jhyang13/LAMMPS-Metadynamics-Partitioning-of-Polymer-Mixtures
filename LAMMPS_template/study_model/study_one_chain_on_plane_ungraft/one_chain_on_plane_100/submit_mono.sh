#!/bin/bash 

NP=1 # the number of processors

	for N in 100 # number of monomers per chain 
	do 
		for EPS in 0.95 1.00 1.05 1.15 
		do
			
			DIR="N${N}_EPS${EPS}"
			mkdir $DIR
			cp templates/def.chain templates/poly.data templates/in.mono templates/submit.mino $DIR
			cd $DIR
                                                    
			sed '{s/{EPS}/'"$EPS"'/g}' < in.mono > tmp1
			mv tmp1 in.run
			rm tmp* in.mono

			sed '{s/{NP}/'"$NP"'/g}' < submit.mino > tmp1
			sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
			mv tmp2 submit.mino
			rm tmp*	
											
			qsub submit.mino
			cd ..

		done
	done
