# Makefile for serial tools

#
# Targets
#

all:
	$(MAKE) chain_orig #micelle2d

#binary2txt:	binary2txt.o
#	g++ -g binary2txt.o -o binary2txt

chain_orig:	chain_orig.o
	ifort chain_orig.o -o chain_orig

#micelle2d:	micelle2d.o
#	ifort micelle2d.o -o micelle2d

#thermo_extract:	thermo_extract.o
#	gcc -g thermo_extract.o -o thermo_extract

clean:
	rm chain_orig #micelle2d
	#rm thermo_extract
	rm *.o

#
# Rules
#

.cpp.o:
	g++ -g -c $<

.c.o:
	gcc -g -c $<

.f.o:	
	ifort -O -w -c $<
