#!/bin/bash 

NP=9 # the number of processors on a node requested to run the jobs submitted here
NSETS=1 # the number of different "species" of polymer chains; 1 = monodisperse, all chains length NA
RUN=1000 # the length of the run in number of timesteps  
FRACA=1 # mol fraction of species 

for X in 5 # x length of box
do
	for Y in 5 # y length of box 
	do
		for Z in 10 # z length of box
		do
			for NA in 1 # number of monomers per chain 
			do
				for PHI in 0.2 # number of monomers per volume
				do 
					for REP in $(seq -f "%03g" 1 1) # an id which distinguishes different replicates of the same run 
					do
						for DATE in 20211213 # the date of when you're submitting the run
						do
									
										DIR="${DATE}_X${X}_Z${Z}_NA${NA}_PHI${PHI}_${REP}"
										mkdir $DIR
										cp templates/def.chain templates/in.lj templates/submit.mino $DIR
										cd $DIR

										RAND=$RANDOM

										sed '{s/{PHI}/'"$PHI"'/g}' < def.chain > tmp1
										sed '{s/{LENGTH}/'"$X"'/g}' < tmp1 > tmp2  
										sed '{s/{HEIGHT}/'"$Z"'/g}' < tmp2 > tmp3  
										sed '{s/{RAND}/'"$RAND"'/g}' < tmp3 > tmp4  
										sed '{s/{NSETS}/'"$NSETS"'/g}' < tmp4 > tmp5  
										sed '{s/{FRACA}/'"$FRACA"'/g}' < tmp5 > tmp6  
										sed '{s/{NA}/'"$NA"'/g}' < tmp6 > tmp7
										mv tmp7 def.chain  
										./../templates/chain_orig < def.chain > poly.data
										rm tmp*

										sed '{s/{HBOX}/'"$X"'/g}' < in.lj > tmp1
									        sed '{s/{HBOY}/'"$Y"'/g}' < tmp1 > tmp2
                                                                                sed '{s/{HBOZ}/'"$Z"'/g}' < tmp2 > tmp3
									        sed '{s/{RUN}/'"$RUN"'/g}' < tmp3 > tmp4
										mv tmp4 in.run
										rm tmp* in.lj

										sed '{s/{NP}/'"$NP"'/g}' < submit.mino > tmp1
										sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
										mv tmp2 submit.mino
										rm tmp*	
											
										qsub submit.mino
										cd ..

	 						done			
						done
					done
				done
			done
		done
	done
