
Files/programs associated with making input coordinates, bonds, etc. of polymer chains:


1.The First Input Files:
   
   def.chain: input file for program "chain_orig" which contains instructions for parameters like the density of monomers, chainlength, etc. 

   chain_orig.f: source code (in Fortran) of program which makes input datafile of polymer chains
   

2. Makefile: makefile which includes instructions for compiliing "chain_orig.f" into binary executable "chain_orig", language translator
  
   #command: make


3. Output files:

   poly.data: output file created by program "chain_orig" combining with "def.chain" which contains atomic coordinates, bonds, box dimensions of polymer chains

   #command: ./chain_orig<def.chain>poly.data




