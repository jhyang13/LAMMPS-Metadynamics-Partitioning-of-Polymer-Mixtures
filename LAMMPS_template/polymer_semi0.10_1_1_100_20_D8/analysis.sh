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

                                                                                	NACHAINS=$(awk 'NR==1{print $9}' $DIR/poly.data)
                                                                                    	NBCHAINS=$(awk 'NR==2{print $9}' $DIR/poly.data)

                                                                                      	cd $DIR

                                                                                      	rm -rf long/
                                                                                        rm -rf short/

                                                                                        mkdir long
                                                                                        mkdir short

  											cd long/
                                                                                    	cp ../../cfg_long.txt .

                                                                                      	sed '{s/{DIR}/'"$DIR"'/g}' < cfg_long.txt > tmp1
                                                                                      	sed '{s/{NA}/'"$NA"'/g}' < tmp1 > tmp2
                                                                                        sed '{s/{NACHAINS}/'"$NACHAINS"'/g}' < tmp2 > tmp3
                                                                                     	sed '{s/{PHIA}/'"$PHIA"'/g}' < tmp3 > tmp4
                                                                                        sed '{s/{EPS}/'"$EPS"'/g}' < tmp4 > tmp5
											sed '{s/{R}/'"$R"'/g}' < tmp5 > tmp6

    											mv tmp6 cfg_long.txt
                                                                                     	rm tmp*
                                                                                 	./../../matrix70 cfg_long.txt                                                                                    			    rm cfg_long.txt
							
											cd ..

                                                                      			cd short/
                                                                                        cp ../../cfg_short.txt .

											sed '{s/{DIR}/'"$DIR"'/g}' < cfg_short.txt > tmp1
                                                                                     	sed '{s/{NB}/'"$NB"'/g}' < tmp1 > tmp2
                                                                                    	sed '{s/{NBCHAINS}/'"$NBCHAINS"'/g}' < tmp2 > tmp3
                                                                                   	sed '{s/{PHIB}/'"$PHIB"'/g}' < tmp3 > tmp4
                                                                                    	sed '{s/{EPS}/'"$EPS"'/g}' < tmp4 > tmp5
    											sed '{s/{R}/'"$R"'/g}' < tmp5 > tmp6                
		
                                                                                        mv tmp6 cfg_short.txt
                                                                                       	rm tmp*
                                                                                   	./../../matrix70 cfg_short.txt
                                                                                   	rm cfg_short.txt
							
											cd ..

											touch ANALYZED

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
