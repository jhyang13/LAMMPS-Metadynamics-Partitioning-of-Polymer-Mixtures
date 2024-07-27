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
				for R in 4
				do
					for D in 8
					do
						for EPS in 0.40 0.60 0.80
						do 
							for DATE in 20230418
							do
								for REP in 01 02 03 04 05 06 07 08 09 10
								do
									
									DIR="${DATE}_NA${NA}_PHIA${PHIA}_NB${NB}_PHIB${PHIB}_R${R}_EPS${EPS}_REP${REP}"

									if [ -d "$DIR" ]
                                                                        then
                                                                        	if [ -e "$DIR/COMPLETED" ]
                                                                                then
                                                                                	if [ $REP != 01 ]
                                                                                        then
												sed '1,1d' $DIR/long/long.density > $DIR/long/long.density.rmvd_first_line
                                                                                                cat $DIR/long/long.density.rmvd_first_line >> PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.long

                                                                                         	sed '1,1d' $DIR/short/short.density > $DIR/short/short.density.rmvd_first_line
                                                                                         	cat $DIR/short/short.density.rmvd_first_line >> PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.short
                                                                                 	else
                                                                                              	cat $DIR/long/long.density >> PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.long
                                                                                                cat $DIR/short/short.density >> PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.short

											fi
										fi
									fi
								done
	 						done			
						done
					done
				done
			done
		done
	done
done
