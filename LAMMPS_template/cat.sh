#!/bin/bash 

for NA in 120
do
	for PHIA in 0.05
	do
		for NB in 40
		do	
			for PHIB in 0.05
			do
				for R in 3
				do
					for D in 6
					do
						for EPS in 0.20 0.40 0.60 0.80 1.00 1.20 1.40 1.60 1.80 2.00
						do 
							for DATE in 20221119
							do
								for REP in 01 02 03
								do
									
									#cat long_contacts_in_pore_EPS${EPS}_${REP}.txt > long_contacts_120_40_1_3_EPS${EPS}.txt
									cat short_contacts_in_pore_EPS${EPS}_${REP}.txt > short_contacts_120_40_1_3_EPS${EPS}.txt 

								done
	 						done			
						done
					done
				done
			done
		done
	done
done
