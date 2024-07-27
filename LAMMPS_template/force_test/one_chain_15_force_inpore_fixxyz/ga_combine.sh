#!/bin/bash 

NP=1

for EPS in 1.00
do
	touch comz_potential_EPS${EPS}_lad10.txt

	for N in 15
	do	
		for R in 2
		do
			for D in 4
			do
				for DATE in 20230426
				do
					for h in -20 -15 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0
					do
									
						DIR="${DATE}_N${N}_R${R}_EPS${EPS}_h${h}_lad10"
						cd $DIR

						cat ana.dat >> /home/jhyang/single_chain_15_force_inpore_fixxyz/comz_potential_EPS${EPS}_lad10.txt
						cd ..

					done
	 			done			
			done
		done
	done
done
