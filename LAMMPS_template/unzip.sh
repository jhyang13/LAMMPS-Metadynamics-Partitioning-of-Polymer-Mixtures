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
						for EPS in 1.80 2.00
						do 
							for DATE in 20221119
							do
								for REP in 01 02 03
								do									
									DIR="${DATE}_NA${NA}_PHIA${PHIA}_NB${NB}_PHIB${PHIB}_R${R}_EPS${EPS}_REP${REP}"

									if [ -d "$DIR" ]
                                                                   	then
                                                                        	if [ -e "$DIR/COMPLETED" ]
                                                                           	then


                                                                                      	cd $DIR

                                                                                        gzip -d long.lammpstrj.gz
											gzip -d short.lammpstrj.gz

                                                                                    	cd ..

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
