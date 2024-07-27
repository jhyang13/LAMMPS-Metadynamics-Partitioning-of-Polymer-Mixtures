#!/bin/bash 

NP=9 # the number of processors on a node requested to run the jobs submitted here
NSETS=2 # the number of different "species" of polymer chains; 1 = monodisperse, all chains length NA
RUN=10000000 # the length of the run in number of timesteps  

for L in 15 # the length of the edges of the simulation box in the x and y dimensions
do
	for H in 15 # the height, or length of the edge, of the simulation box in the z dimension
	do
		for NA in 6 # the number of monomers per polymer chain 
		do
			for PHIA in 0.35 # the monomer number concentration of the polymer chains, i.e., #monomers/(volume above pore)
			do
				for R in 1.0 1.2 1.4 # the radius of the pore; note that the effective radius is less than this because of potential used
				do
					for D in 8 # the length of the pore, from opening to bottom 
					do
						for EPS in 1.0 # the strength of the interaction between a polymer monomer and 
						do 
							for REP in $(seq -f "%03g" 1 1) # an id which distinguishes different replicates of the same run 
							do
								for DATE in 20211004 # the date of when you're submitting the run; makes organizing many runs easy
								do
									
											DIR="${DATE}_L${L}_H${H}_NA${NA}_PHIA${PHIA}_R${R}_D${D}_EPS${EPS}_${REP}" # the name of the folder created; 1 job = 1 folder 
											mkdir $DIR
											cp templates/def.chain templates/in.mono templates/submit.quest $DIR
											cd $DIR

											PHI=$PHIA
											FRACA=1
											RAND=$RANDOM

											sed '{s/{PHI}/'"$PHI"'/g}' < def.chain > tmp1
											sed '{s/{LENGTH}/'"$L"'/g}' < tmp1 > tmp2  
											sed '{s/{HEIGHT}/'"$H"'/g}' < tmp2 > tmp3  
											sed '{s/{RAND}/'"$RAND"'/g}' < tmp3 > tmp4  
											sed '{s/{NSETS}/'"$NSETS"'/g}' < tmp4 > tmp5  
											sed '{s/{FRACA}/'"$FRACA"'/g}' < tmp5 > tmp6  
											sed '{s/{NA}/'"$NA"'/g}' < tmp6 > tmp7
											mv tmp7 def.chain  
											./../templates/chain_orig < def.chain > poly.data
											rm tmp*
											
											HBOX=$(bc -l <<< "$L/2")
											ZLO=$(bc -l <<< "$D+1")
											ZHI=$(bc -l <<< "$H")
											RADHI=$(bc -l <<< "$R+$R")
											LO=$(bc -l <<< "$RADHI-$R")
											RAND1=$RANDOM
											RAND2=$RANDOM
											RAND3=$RANDOM
											RAND4=$RANDOM
											sed '{s/{RAND1}/'"$RAND1"'/g}' < in.mono > tmp1
											sed '{s/{HBOX}/'"$HBOX"'/g}' < tmp1 > tmp2
											sed '{s/{EPS}/'"$EPS"'/g}' < tmp2 > tmp3
											sed '{s/{D}/'"$D"'/g}' < tmp3 > tmp4
											sed '{s/{R}/'"$R"'/g}' < tmp4 > tmp5
											sed '{s/{RADHI}/'"$RADHI"'/g}' < tmp5 > tmp6
											sed '{s/{ZLO}/'"$ZLO"'/g}' < tmp6 > tmp7
											sed '{s/{ZHI}/'"$ZHI"'/g}' < tmp7 > tmp8
											sed '{s/{LO}/'"$LO"'/g}' < tmp8 > tmp9
											sed '{s/{RAND2}/'"$RAND2"'/g}' < tmp9 > tmp10
											sed '{s/{RAND3}/'"$RAND3"'/g}' < tmp10 > tmp11
											sed '{s/{RAND4}/'"$RAND4"'/g}' < tmp11 > tmp12
											sed '{s/{RUN}/'"$RUN"'/g}' < tmp12 > tmp13
											mv tmp13 in.run
											rm tmp* in.mono

											sed '{s/{NP}/'"$NP"'/g}' < submit.quest > tmp1
											sed '{s/{DIR}/'"$DIR"'/g}' < tmp1 > tmp2
											mv tmp2 submit.quest
											rm tmp*	
											
											sbatch submit.quest 
											cd ..

	 								done			
								done
							done
						done
					done
				done
			done
		done
	done
