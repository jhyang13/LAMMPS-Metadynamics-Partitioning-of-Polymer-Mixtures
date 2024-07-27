#!/bin/bash 

rm long_partition.txt
rm short_partition.txt

touch long_partition.txt
touch short_partition.txt

for PHIA in 1
do	
	for PHIB in 1
	do
		for EPS in 0.40 0.60 0.80
		do 	
			 rm ana.dat
			 ga PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.long
          		 cat ana.dat >> /home/jhyang/polymer_semi0.10_1_1_100_20_D8/long_partition.txt

                         rm ana.dat                         
			 ga PHIA${PHIA}_PHIB${PHIB}_EPS${EPS}.partcoeff.short
                         cat ana.dat >> /home/jhyang/polymer_semi0.10_1_1_100_20_D8/short_partition.txt

		done
	done
done
