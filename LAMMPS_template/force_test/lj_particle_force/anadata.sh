#!/bin/bash 

NP=1
EPS=1.00

touch force.txt
touch etotal.txt

for N in 1
do	 
	for DATE in 20230406
	do
		for h in 0.80 0.85837 0.90 0.92 0.94 0.96 0.98 1.00 1.20 1.40 1.60 1.80 2.00 2.20 2.40 2.60
		do 

			DIR="${DATE}_N${N}_EPS${EPS}_h${h}_ljforce_lan10"
			cd $DIR
			
			rm ana.dat							
			ga comz_forcez_EPS${EPS}.txt
                        cat ana.dat >> /home/jhyang/lj_particle_force/force.txt
                        
			rm ana.dat
			ga comz_etotal_EPS${EPS}.txt
			cat ana.dat >> /home/jhyang/lj_particle_force/etotal.txt

			cd ..

		done
	done
done
