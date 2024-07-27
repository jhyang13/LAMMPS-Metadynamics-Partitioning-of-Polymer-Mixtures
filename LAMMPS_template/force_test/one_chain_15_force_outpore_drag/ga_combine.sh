#!/bin/bash 

NP=1

for EPS in 3.00 4.00 5.00
do
	touch comz_forcez_EPS${EPS}_lad100.txt

	for N in 15
	do	
		for R in 2
		do
			for D in 4
			do
				for DATE in 20230328
				do
					for h in 1 2 3 4 5 6 7 8 9 10 15 20
					do
									
						DIR="${DATE}_N${N}_R${R}_EPS${EPS}_h${h}_lad100"
						cd $DIR

						cat ana.dat >> /home/jhyang/single_chain_15_force_outpore_fixxyz/comz_forcez_EPS${EPS}_lad100.txt
						cd ..

					done
	 			done			
			done
		done
	done
done
