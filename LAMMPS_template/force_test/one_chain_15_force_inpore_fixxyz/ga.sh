#!/bin/bash 

NP=1

for N in 15
do	
	for R in 2
	do
		for D in 4
		do
			for EPS in 1.00
			do 
				for DATE in 20230426
				do
					for h in -20 -15 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0
					do
									
						DIR="${DATE}_N${N}_R${R}_EPS${EPS}_h${h}_lad10"
						cd $DIR
						rm ana.dat
						ga comz_potential_EPS${EPS}.txt
						cd ..

					done
	 			done			
			done
		done
	done
done
