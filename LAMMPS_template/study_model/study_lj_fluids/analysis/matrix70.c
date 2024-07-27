/*
  Dump analysis for LAMMPS 2004
  Won-Ki Roh, Dec 2004/Jan 2005
  Modifications EL
  * revised Feb 1 2005 (add function to calculate bond length) *
  * revised Feb 4 2005 (add function to measure no. of adsorbed monomers)*v.54
  * revised Mar 11 2005 (add function to measure no. of adsorbed chains)*v.56
  * revised Mar 11-12 2005 (add function to measure no. of loop,train,tail's monomers)*v.57
  * revised Mar 22-24 2005 (add function to measure length & no. of loop,train,tail)*v.58
  * revised Mar 26    2005 (add function to calculate error bar of msd)*v.59
  * revised Apr 8     2005 (checking the dump file's disorder of configuration)*v.59_2
  * revised Apr 13    2005 (fix error bar of msd)*v.59_3
  Modifications Max Meirow
  * revised Aug 19    2020 (update to read LAMMPS 12 Dec 2018 dump files)*v.60
  * revised Aug 26    2020 (update to calculate number density and vol. fraction in bulk and pore)*v.61
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc2d.h"
#include "alloc3d.h"

#define TRUE  1
#define FALSE 0


#define KEYWORDS {"dumpname", "firsttimestep", "saveevery", "timesteps", "maxdeltat", "rho", "chainlength", "no_chains", "z_com", "length", "height", "depth", "radius", "bin_length", "bin_width", "bin_cutoff", "disk_radius", "disk_height", "ad_condition_z", "compute_msd", "outputfile", "msdout", "adsorb", "aveforce", "centerofmass", "ejectiontime", "pore", "gofz", "densities", "energies", "opening", "bottom_disk_height", "bottom", "partition"}
#define N_KEYWORDS 34

FILE *fp;
FILE *out1;
FILE *out2;
FILE *out3;
FILE *out4;
FILE *out5;
FILE *out6;
FILE *out7;
FILE *out8;
FILE *out9;
FILE *out10;
FILE *out11;
FILE *out12;
FILE *out13;

int main(int argc, char *argv[])
{
	FILE *cfg;
	char keyword[1024], value[1024], dumpname[1024], command[1024], outputfile[1024], msdout[1024], adsorb[1024], aveforce[1024], centerofmass[1024], ejectiontime[1024], pore[1024], gofz[1024], densities[1024], energies[1024], opening[1024], bottom[1024], partition[1024];
	char *keywordlist[] = KEYWORDS;
	double **position_array, ***cm_array;
	int timestep, firsttimestep, saveevery, i, k, checked, lineno, readstatus, compressed, length;
	int notimesteps, maxdelta, chainlength, nochains, msdon, inpore;	
	int total_nomon_ad, total_nochain_ad, no_train_mon, no_tail_mon, no_loop_mon, no_nonadsorb_mon, no_train, no_tail, no_loop, no_nonadsorb;
	int atomno, monomer, chain, monomers[2];
	double avetrain_mon, avetail_mon, aveloop_mon, avenonadsorb_mon;
	double ad_condition=0.0;
	double zcom, rho;
	double aveRg, aveaveRg, aveRgsum, aveRe, aveavel, coordaveRg[3], coordaveRe[3], coordaveavel[3];
	double fave, favecomp[3], favesum, avefave;
	double numdensity[5], partcoeff;
	double edgelength, height, depth, radius;
	double binlength, binwidth, bincut, center;
	int nbins, nxbins, nybins, j;
	double dist, xval, yval;
	double diskrad, diskheight, bdiskheight;
	
	
	int readconfig(double **positions, int nchains, int chainlength);
	void printconfig(double **positions, int atoms);
	void calcbottom(double **positions, int chainlength, int nchains, double depth, double bdiskheight, double rho);
	void calcopening(double **positions, int chainlength, int nchains, double edgelength, double height, double diskrad, double diskheight, int firsttimestep, int timestep, int saveevery);
	void calcbindensity(double **positions, int chainlength, int nchains, double edgelength, double binlength);
	void calczdist(double **positions, int chainlength, int nchains, double edgelength, double height, int firsttimestep, int timestep, int saveevery, double binwidth, double bincut);
	void calcraddist(double **positions, int chainlength, int nchains, double radius, int firsttimestep, int timestep, int saveevery);
	void chaininpore(double **positions, int chainlength, int nchains, int firsttimestep, int timestep, int saveevery);
	double calcforce(double **positions, int chainlength, int nchains, double favecomp[3]);
	double calcgyration(double ***cm_array, double **positions, int chainlength, int nchains, int timestep, double coordaveRg[3]);
	void calccentermass(double ***cm_array, double **positions, int chainlength, int nchains, int timestep);
	double calcmonomerdistance(double **positions, int chainlength, int nchains, double aveavel, double coordaveavel[3]);
	double calcnumdensity(double **positions, int chainlength, int nchains, double edgelength, double height, double depth, double radius, double partcoeff, double numdensity[5], int monomers[2]);
	double calcendlength(double **positions, int chainlength, int nchains, double aveRe, double coordaveRe[3]);
	void calcmsd(double ***cm_array, int notimesteps, int nchains, int maxdelta);
	int adsorptioncondition(double **positions, int chainlength, int nchains, int total_nomon_ad, double ad_condition, int *total_nochain_ad);
	void looptraintail(double **positions, int chainlength, int nchains, double ad_condition, int *no_train_mon, int *no_tail_mon, int *no_loop_mon, int *no_nonadsorb_mon, double *avetrain_mon, double *avetail_mon, double *aveloop_mon, double *avenonadsorb_mon, int *no_train, int *no_tail, int *no_loop, int *no_nonadsorb);
	double calcatomdist(double **positions);

	if (argc<2)
	{
		fprintf(stderr, "Usage: %s configurationfile\n", argv[0]);
		exit(0);
	}
	
	cfg=fopen(argv[1], "r");

	if (cfg == NULL)
	{
		fprintf(stderr, "Cannot open configuration file %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stderr, "Reading configuration file %s\n", argv[1]);
	}

	readstatus = fscanf(cfg, "%s = %s\n", keyword, value);
	lineno = 1;

	while (readstatus != EOF)
	{
	    i=0;
	    while (i<N_KEYWORDS)
	    {
		if (strcmp(keyword, keywordlist[i])==0) break;
		i++;
	    }
	    if (i == N_KEYWORDS)
	    {
		fprintf(stderr, "Parse error in line %d\n", lineno);
		exit(EXIT_FAILURE);
	    }
	    switch(i)
	    {
		case(0):
		    strcpy(dumpname, value);
		    printf("\n");
		    printf("Dumpname: %s\n", dumpname);
		    break;
		case(1):
		    firsttimestep=atoi(value);
		    printf("First timestep: %d\n", firsttimestep);
		    break; 
		case (2):
		    saveevery=atoi(value);
		    printf("Save every n timesteps: %d\n", saveevery);
		    break;    
		case(3):
		    notimesteps=atoi(value);
		    printf("No. of timesteps: %d\n", notimesteps);
		    break;
		case(4):
		    maxdelta=atoi(value);
		    printf("Max. delta t: %d\n", maxdelta);
		    break;
		case(5):
		    rho=atof(value);
		    printf("Monomer number density: %.2f\n", rho);
		    break;    
		case(6):
		    chainlength=atoi(value);
		    printf("Chainlength: %d\n", chainlength);
		    break;
		case(7):
		    nochains=atoi(value);
		    printf("No. of chains: %d\n", nochains);
		    break;
		case(8):
		    zcom=atof(value);
		    printf("Z-coordinate center-of-mass: %10.5f\n", zcom);
		    break;    
		case(9):
		    edgelength=atof(value);
		    printf("Bulk length: %10.5f\n", edgelength);
	       	    break;
	        case(10):
	            height=atof(value);
	            printf("Bulk height: %10.5f\n", height);
		    break;
	        case(11):
		    depth=atof(value);
		    printf("Pore depth: %10.5f\n", depth);
		    break;
	        case(12):
		    radius=atof(value);
		    printf("Pore radius: %10.5f\n", radius);
		    break;	
		case(13):
		    binlength=atof(value);
		    printf("Length of cubic bins: %10.5f\n", binlength);
		    break;    
		case(14):
		    binwidth=atof(value);
	            printf("Width of g(z) bin: %10.5f\n", binwidth);
		    break;
		case(15):
		    bincut=atof(value);
		    printf("Upper cutoff for g(z): %10.5f\n", bincut);
		    break;	  
		case(16):
		    diskrad=atof(value);
		    printf("Radius of disk above opening: %10.5f\n", diskrad);
		    break;
		case(17):
		    diskheight=atof(value);
		    printf("Height of disk above opening: %10.5f\n", diskheight);
		    break;    
		case(18):
		    ad_condition=atof(value);
		    printf("Adsorb condition z: %10.5f\n", ad_condition);
		    break;
		case(19):
		    msdon=atoi(value);
		    if(msdon==1)
			printf("Compute msd.\n");
		    else
			printf("Don't compute msd.\n");
		    break;
		case(20):
		    strcpy(outputfile, value);
		    printf("outputfile -> %s\n", outputfile);
		    break;
		case(21):
		    strcpy(msdout, value);
		    printf("msdout -> %s\n", msdout);
		    break;
		case(22):
		    strcpy(adsorb, value);
		    printf("adsorb -> %s\n", adsorb);
		    break;
		case(23):
		    strcpy(aveforce, value);
		    printf("aveforce -> %s\n", aveforce);
		    break;    
		case(24):
		    strcpy(centerofmass, value);
		    printf("centerofmass -> %s\n", centerofmass);
		    break;
		case(25):
		    strcpy(ejectiontime, value);
		    printf("ejectiontime -> %s\n", ejectiontime);
		    break;    
		case(26):
		    strcpy(pore, value);
		    printf("pore data -> %s\n", pore);
		    break;
		case(27):
		    strcpy(gofz, value);
		    printf("g(z) vs. time data -> %s\n", gofz);
		    break;
		case(28):
		    strcpy(densities, value);
		    printf("bin densities vs. time -> %s\n", densities);
		    break;
		case(29):
		    strcpy(energies, value);
		    printf("bin energies vs. time -> %s\n", energies);
		    break;
		case(30):
		    strcpy(opening, value);
		    printf("Monomers near pore opening -> %s\n", opening);
		    break;
		case(31):
		    bdiskheight=atof(value);
		    printf("Height of disk at bottom of pore -> %10.5f\n", bdiskheight);
		    break;
		case(32):
		    strcpy(bottom, value);
		    printf("Chain ends near bottom -> %s\n", bottom);
		    break;    
		case(33):
				strcpy(partition, value);
				printf("Partitioning output file -> %s\n", partition);
				break;

		default:
		    fprintf(stderr, "Internal error\n");
		    exit(EXIT_FAILURE);
	    }			
	    readstatus = fscanf(cfg, "%s = %s\n", keyword, value);
	    lineno++;
	}
	printf("\n");
	
	fclose(cfg);
	
	position_array = allocate_double_matrix(nochains*chainlength, 7);
	cm_array = allocate_3d_double_matrix(notimesteps,  nochains, 3);
	
        length = strlen(dumpname);
        if (length > 3)
        {
                if (strcmp(&dumpname[length-3], ".gz") == 0)
                        compressed=TRUE;
                else
                        compressed=FALSE;
        }
        else
        {
                compressed=FALSE;
        }
	if ( compressed == TRUE )
        {
                strcpy(command, "gzip -dc ");
                strcat(command, dumpname);
                fp = popen(command, "r");
        }
        else
        {
                fp  = fopen(dumpname, "r");
        }

	if (fp==NULL)
	{
		fprintf(stderr, "Cannot open dumpfile %s\n", dumpname);
		exit(EXIT_FAILURE);
	}

	
	out1=fopen(outputfile, "w");
  out13=fopen(partition, "w");
	//out8=fopen(gofz, "w");
	//out9=fopen(densities, "w");
	//out10=fopen(energies, "w");
	//out11=fopen(opening, "a");
	//out12=fopen(bottom, "a");

	/*out4=fopen(aveforce, "a");*/

	/*out5=fopen(centerofmass, "w");
	out6=fopen(ejectiontime, "a");
	out7=fopen(pore, "w");*/
	
		
	fprintf(out1, "# timestep  xRg^2        yRg^2           zRg^2           Rg^2          xRe^2         yRe^2         zRe^2          Re^2          xL^2           yL^2          zL^2           L^2\n");
	
	fprintf(out13, "# timestep phi_pore phi_bulk part_coeff\n");

	//fprintf(out11, "# timestep monomer_count monomer_density norm_density end_count end_density norm_enddensity\n");
	//fprintf(out1, "# timestep phi_pore phi_bulk part_coeff\n");
	
	/*fprintf(out5, "# timestep   xCOM	 yCOM		 zCOM\n");
	fprintf(out7, "# timestep  monomersinpore\n");*/

/*	out3=fopen(adsorb, "w");
	fprintf(out3, "# tstep admon(/%d) ad_percent adchains trainmon tailmon loopmon nonadmon ave_trainmon ave_tailmon ave_loopmon ave_nonadmon no_train  no_tail   no_loop\n", chainlength*nochains);
*/
	aveRgsum=0;
	favesum=0;
	checked=0;

	/*// write the centers of the histogram bins 
	
	nbins = (int) bincut/binwidth;
	double bincenters[nbins];
	fprintf(out8, "# ");
	for(j=0; j < nbins; j++)
	{
		center = j*binwidth + binwidth/2.0;
		bincenters[j] = center;
		fprintf(out8, "%10.4f ", bincenters[j]);
	}	
	fprintf(out8, "\n");

	nxbins = (int) edgelength/binlength;
	nybins = (int) edgelength/binlength;
	fprintf(out9, "# ");
	fprintf(out10, "# ");
	for(i = 0; i < nxbins; i++)
	{
		for(j = 0; j < nybins; j++)
		{
			xval = (-1*(edgelength/2.0)) + (binlength/2.0) + (i*binlength);
			yval = (-1*(edgelength/2.0)) + (binlength/2.0) + (j*binlength);
			
			fprintf(out9, "(%3.1f, %3.1f) ", xval, yval);
			fprintf(out10, "(%3.1f, %3.1f) ", xval, yval);
		}
	}
	fprintf(out9, "\n");
	fprintf(out10, "\n");
	*/

	for(timestep=0; timestep<notimesteps; timestep++)
	{
	    
	    readconfig(position_array, nochains, chainlength);
	    
	    /*printconfig(position_array, nochains*chainlength);*/
	    
	    calccentermass(cm_array, position_array, chainlength, nochains, timestep);

	    aveavel=calcmonomerdistance(position_array, chainlength, nochains, aveavel, coordaveavel);

	    aveRe=calcendlength(position_array, chainlength, nochains, aveRe, coordaveRe);	
	    
	    aveRg=calcgyration(cm_array, position_array, chainlength, nochains, timestep, coordaveRg);
	   
    	    //calczdist(position_array, chainlength,nochains, edgelength, height, firsttimestep, timestep, saveevery, binwidth, bincut);
	    //calcbindensity(position_array, chainlength, nochains, edgelength, binlength);
	    //dist = calcatomdist(position_array);
	    
	    /*if(timestep == 0)
	    {
		calcopening(position_array, chainlength, nochains, edgelength, height, diskrad, diskheight, firsttimestep, timestep, saveevery);
	    }

	    if(timestep == notimesteps - 1)
	    {
		calcbottom(position_array, chainlength, nochains, depth, bdiskheight, rho);	
	    }*/
	    
	    /*if(timestep % 100 == 0){

	    	calcraddist(position_array, chainlength, nochains, radius, firsttimestep, timestep, saveevery);
	    }*/
	    //chaininpore(position_array, chainlength, nochains, firsttimestep, timestep, saveevery);

	    /*fave=calcforce(position_array, chainlength, nochains, favecomp);*/	

	    partcoeff=calcnumdensity(position_array, chainlength, nochains, edgelength, height, depth, radius, partcoeff, numdensity, monomers);
	    
	    /*looptraintail(position_array, chainlength, nochains, ad_condition, &no_train_mon, &no_tail_mon, &no_loop_mon, &no_nonadsorb_mon, &avetrain_mon, &avetail_mon, &aveloop_mon, &avenonadsorb_mon, &no_train, &no_tail, &no_loop, &no_nonadsorb);*/

	    /*printf("%d %d %d %d\n", nonadsorb, train, tail, loop);*/

	    /*total_nomon_ad=adsorptioncondition(position_array, chainlength, nochains, total_nomon_ad, ad_condition, &total_nochain_ad);
	 
	    fprintf(out3, "    %d     %d %14.4f      %d        %d       %d        %d       %d    %10.4f  %10.4f  %10.4f   %10.4f %10.4f %10.4f %10.4f\n", timestep, total_nomon_ad, (double)total_nomon_ad/(chainlength*nochains), total_nochain_ad, no_train_mon, no_tail_mon, no_loop_mon, no_nonadsorb_mon, avetrain_mon, avetail_mon, aveloop_mon, avenonadsorb_mon, (double)no_train/nochains, (double)no_tail/nochains, (double)no_loop/nochains);*/

	    fprintf(out13, "  %d %14.10f %14.10f %14.10f\n", firsttimestep+timestep*saveevery, numdensity[0], numdensity[1], partcoeff);

	    fprintf(out1, "  %d %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", firsttimestep+timestep*saveevery, coordaveRg[0], coordaveRg[1], coordaveRg[2], aveRg, coordaveRe[0], coordaveRe[1], coordaveRe[2], aveRe, coordaveavel[0], coordaveavel[1], coordaveavel[2], aveavel);

	    /*fprintf(out5, "  %d %14.10f %14.10f %14.10f\n", timestep*saveevery, cm_array[timestep][0][0], cm_array[timestep][0][1], cm_array[timestep][0][2]);

	    fprintf(out7, "  %d %d\n", firsttimestep+timestep*saveevery, monomers[0]);*/ 
	    
	    aveRgsum+=aveRg;
	    
	   /* if(timestep > 999)
	    {
		favesum+=fave;	
	    }

	    if(checked > 0)
	    {
		    continue;
	    }	
	    else 
	    {
            	k=0;
		for(i=0; i<chainlength; i++)
		{
			if(position_array[i][2] < 0)
			{
				k++;
			}
		}

		if(k < 1)
		{
			fprintf(out6, "  %d %d\n", chainlength, timestep*saveevery);
			checked++;
		}
	    }*/	   
	    
	}

	


	aveaveRg=aveRgsum/notimesteps;
	/*printf("%d\n", notimesteps-1000);*/
	/*avefave=favesum/(notimesteps-1000);*/
	/*printf("Average force over last timesteps: %14.10f\n", avefave);*/   
	
	/*fprintf(out4, "%lg %lg\n", zcom, avefave);*/

       	if(msdon==1)
	{
	    out2=fopen(msdout, "w");
	    fprintf(out2, "# delta t    msd(x)        error            msd(y)      error            msd(z)         error             mad        error          msd(x+y)      error         no.timesteps(%d)\n", notimesteps);
	    
	    calcmsd(cm_array, notimesteps, nochains, maxdelta);
	    fclose(out2);
	}
		
	if (compressed == TRUE)
	{
	    pclose(fp);
	}
	else
	{
	    fclose(fp);
	}
	
	fclose(out1);
	fclose(out13);
	//fclose(out8);
	//fclose(out9);
	//fclose(out10);
	//fclose(out11);
	//fclose(out12);
	/*fclose(out5);
	fclose(out6);
	fclose(out7);
	fclose(out3);*/

	free_3d_double_matrix(cm_array);
	free_double_matrix(position_array);

	exit(0);
}

void calcmsd(double ***cm_array, int notimesteps, int nchains, int maxdelta)
{
	int i, j, delta, chain, deltat, nodeltas, max_j, max_nodeltas; 
	double disp[3], **sumd, **msd, **xymsd, ***coordmsd; 
	double **sqrsumd, ***sqrcoordmsd, **sqrmsd, **sqrxymsd;
	double summsd, sumxymsd, coordsummsd[3]; 
	double sqrsummsd, sqrsumxymsd, sqrcoordsummsd[3];
	double varmsd, varxymsd, varcoordmsd[3]; 
	double errmsd, errxymsd, errcoordmsd[3];
	max_nodeltas=notimesteps-1;
	double number_of_samples[max_nodeltas];	
	double tempmsd[max_nodeltas], tempxymsd[max_nodeltas], tempsqrmsd[max_nodeltas], tempsqrxymsd[max_nodeltas], nodelta[max_nodeltas];
       

	sumd = allocate_double_matrix(notimesteps, 3);
	sqrsumd = allocate_double_matrix(notimesteps, 3);
	msd = allocate_double_matrix(nchains, notimesteps);
	sqrmsd = allocate_double_matrix(nchains, notimesteps);
	xymsd = allocate_double_matrix(nchains, notimesteps);
	sqrxymsd = allocate_double_matrix(nchains, notimesteps);
	coordmsd = allocate_3d_double_matrix(nchains, notimesteps, 3);
	sqrcoordmsd = allocate_3d_double_matrix(nchains, notimesteps, 3);

	for(chain=0; chain<nchains; chain++)
	{
	    /* Initialize all of the matrix values to zero*/	
	    for(delta=0; delta<maxdelta; delta++)
	    {
		sumd[delta][0]=0;
		sumd[delta][1]=0;
		sumd[delta][2]=0;
		sqrsumd[delta][0]=0;
		sqrsumd[delta][1]=0;
		sqrsumd[delta][2]=0;
		tempmsd[delta]=0;
		tempxymsd[delta]=0;
		tempsqrmsd[delta]=0;
		tempsqrxymsd[delta]=0;
		nodelta[delta]=0;

		number_of_samples[delta]=0;
	    }
	    /*printf("%g\n", number_of_samples[0]);*/

	    for(i=0; i<(notimesteps-1); i++)
	    {
		max_j=i+maxdelta;
		
		if(max_j > (notimesteps-1))
		    max_j=notimesteps-1;
	
		for(j=i; j<max_j; j++)
		{
		    deltat=j-i;
		    disp[0]=cm_array[i][chain][0]-cm_array[j+1][chain][0];
		    disp[1]=cm_array[i][chain][1]-cm_array[j+1][chain][1];
		    disp[2]=cm_array[i][chain][2]-cm_array[j+1][chain][2];
		    
		    sumd[deltat][0]+=disp[0]*disp[0];
		    sumd[deltat][1]+=disp[1]*disp[1];
		    sumd[deltat][2]+=disp[2]*disp[2];

		    sqrsumd[deltat][0]+=disp[0]*disp[0]*disp[0]*disp[0];
		    sqrsumd[deltat][1]+=disp[1]*disp[1]*disp[1]*disp[1];
		    sqrsumd[deltat][2]+=disp[2]*disp[2]*disp[2]*disp[2];
		 
		    tempmsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		    tempxymsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]);

		    tempsqrmsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2])*(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		    tempsqrxymsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1])*(disp[0]*disp[0]+disp[1]*disp[1]);

		    number_of_samples[deltat]++;
		}    
	    }

	    for(delta=0; delta<maxdelta; delta++)
	    { 
		nodeltas=notimesteps-(delta+1);
		
		if(number_of_samples[delta]!=nodeltas)
		 printf("Ooops, It's wrong number of deltat!\n");		
		coordmsd[chain][delta][0]=sumd[delta][0]/nodeltas;
		coordmsd[chain][delta][1]=sumd[delta][1]/nodeltas;
		coordmsd[chain][delta][2]=sumd[delta][2]/nodeltas;

		sqrcoordmsd[chain][delta][0]=sqrsumd[delta][0]/nodeltas;
		sqrcoordmsd[chain][delta][1]=sqrsumd[delta][1]/nodeltas;
		sqrcoordmsd[chain][delta][2]=sqrsumd[delta][2]/nodeltas;
		
		msd[chain][delta]=tempmsd[delta]/nodeltas;
		sqrmsd[chain][delta]=tempsqrmsd[delta]/nodeltas;

		xymsd[chain][delta]=tempxymsd[delta]/nodeltas;
		sqrxymsd[chain][delta]=tempsqrxymsd[delta]/nodeltas;

		nodelta[delta]=nodeltas;
	    } 
	}
	
	for(delta=0; delta<maxdelta; delta++)
	{
	    coordsummsd[0]=0;
	    coordsummsd[1]=0;
	    coordsummsd[2]=0;
	    sqrcoordsummsd[0]=0;
	    sqrcoordsummsd[1]=0;
	    sqrcoordsummsd[2]=0;

	    summsd=0;
	    sumxymsd=0;
	    sqrsummsd=0;
	    sqrsumxymsd=0;
	    
	    for(chain=0; chain<nchains; chain++)
	    {  
		coordsummsd[0]+=coordmsd[chain][delta][0];
		coordsummsd[1]+=coordmsd[chain][delta][1];
		coordsummsd[2]+=coordmsd[chain][delta][2];

		sqrcoordsummsd[0]+=sqrcoordmsd[chain][delta][0];
		sqrcoordsummsd[1]+=sqrcoordmsd[chain][delta][1];
		sqrcoordsummsd[2]+=sqrcoordmsd[chain][delta][2];

		summsd+=msd[chain][delta];
		sumxymsd+=xymsd[chain][delta];

		sqrsummsd+=sqrmsd[chain][delta];
		sqrsumxymsd+=sqrxymsd[chain][delta];

	    }
	    
	    coordsummsd[0]/=nchains;
	    coordsummsd[1]/=nchains;
	    coordsummsd[2]/=nchains;

	    sqrcoordsummsd[0]/=nchains;
	    sqrcoordsummsd[1]/=nchains;
	    sqrcoordsummsd[2]/=nchains;

	    summsd/=nchains;
	    sumxymsd/=nchains;

	    sqrsummsd/=nchains;
	    sqrsumxymsd/=nchains;

	    varcoordmsd[0]=sqrcoordsummsd[0]-coordsummsd[0]*coordsummsd[0];
	    varcoordmsd[1]=sqrcoordsummsd[1]-coordsummsd[1]*coordsummsd[1];
	    varcoordmsd[2]=sqrcoordsummsd[2]-coordsummsd[2]*coordsummsd[2];

	    varmsd=sqrsummsd-summsd*summsd;
	    varxymsd=sqrsumxymsd-sumxymsd*sumxymsd;

	    if (nchains==1)
	    {
	    errcoordmsd[0]=sqrt(varcoordmsd[0])/sqrt(nchains);
	    errcoordmsd[1]=sqrt(varcoordmsd[1])/sqrt(nchains);
	    errcoordmsd[2]=sqrt(varcoordmsd[2])/sqrt(nchains);
	    errmsd=sqrt(varmsd)/sqrt(nchains);
	    errxymsd=sqrt(varxymsd)/sqrt(nchains);
	    }
	    else
	    {
	    errcoordmsd[0]=sqrt(varcoordmsd[0])/sqrt(nchains-1);
	    errcoordmsd[1]=sqrt(varcoordmsd[1])/sqrt(nchains-1);
	    errcoordmsd[2]=sqrt(varcoordmsd[2])/sqrt(nchains-1);
	    errmsd=sqrt(varmsd)/sqrt(nchains-1);
	    errxymsd=sqrt(varxymsd)/sqrt(nchains-1);
	    }
	    
	    /* printf("%g %g %g %g %g\n", sqrcoordsummsd[0], sqrcoordsummsd[1], sqrcoordsummsd[2], sqrsummsd, sqrsumxymsd);
	    printf("%g %g %g %g %g\n", varcoordmsd[0], varcoordmsd[1], varcoordmsd[2], varmsd, varxymsd);
	    printf("  %d  %14.10f %14.10f %14.10f %14.10f %14.10f\n", delta+1, coordsummsd[0], coordsummsd[1], coordsummsd[2], summsd, sumxymsd);
	    printf("%g %g %g %g %g\n", errcoordmsd[0], errcoordmsd[1], errcoordmsd[2], errmsd, errxymsd);*/

	    fprintf(out2, "    %d %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f\n", delta+1, coordsummsd[0], errcoordmsd[0], coordsummsd[1], errcoordmsd[1], coordsummsd[2], errcoordmsd[2], summsd, errmsd, sumxymsd, errxymsd);
	    
	}	
	
	free_double_matrix(sumd);
	free_double_matrix(msd);
	free_double_matrix(xymsd);
	free_3d_double_matrix(coordmsd);
	free_double_matrix(sqrsumd);
	free_double_matrix(sqrmsd);
	free_double_matrix(sqrxymsd);
	free_3d_double_matrix(sqrcoordmsd);

	return;
}

int readconfig(double **positions, int nchains, int chainlength)
{
	int chain, monomer, current_atom;
	int i, headerstatus, timestep, natoms, atomtype;
	double box[3][2], xbox, ybox, zbox;
	
	char line[1024];
	int surface, surfatomno, surfatomtype;
	double surfx, surfy, surfz;

	int readatom(double *atom, double xbox, double ybox, double zbox, int *atomtype);
	int readheader(int *timestep, int *natoms, double box[3][2], int nchains, int chainlength);
	
	headerstatus = readheader(&timestep, &natoms, box, nchains, chainlength);
	if ( headerstatus != 0)
	{
		printf("%d %d %d %d\n", timestep, natoms, nchains, chainlength);
		/*printf("%lg %lg %lg %lg %lg %lg\n", box[0][0], box[0][1], box[1][0], box[1][1], box[2][0], box[2][1]);*/
		printf("Could not read header! Status=%d\n", headerstatus);
		exit(EXIT_FAILURE);
	}
	
	/*printf("\n");    
	printf("Processing timestep: %d\n", timestep);
	printf("Total no of atoms: %d\n", natoms);*/
	xbox=box[0][1]-box[0][0];
	ybox=box[1][1]-box[1][0];
	zbox=box[2][1]-box[2][0];

	
	for(chain=0; chain<nchains; chain++)
	{
		for(monomer=0; monomer<chainlength; monomer++)
		{
			current_atom=chain*chainlength+monomer;
			i = readatom(positions[current_atom], xbox, ybox, zbox, &atomtype);
			/*if (current_atom+1 != i)
			{
			    printf("The configuration of dump is disordered!!!\n Please sort it!!!\n");
			    exit(EXIT_FAILURE);
			}*/
		/*	printf("I want to read atom %d. I actually read atom no. %d\n", current_atom, i);*/
		}
	}
	
	/*fgets(line, 1024, fp);
	puts(line);*/
	
	return(0);
}

int readheader(int *timestep, int *natoms, double box[3][2], int nchains, int chainlength)
{
	/* updated to reflect new headers for BOX BOUND and ATOMS */
	/* Sep 20 2020: updated to include new header for ATOMS */

	int i;
	char line[1024];
	
	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: TIMESTEP\n") != 0)
	{   
                printf("%s\n", line);
		return(-1);
	}
	fscanf(fp, "%d\n", timestep);
	
	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: NUMBER OF ATOMS\n") != 0)
	{
		return(-2);
	}
	fscanf(fp, "%d\n", natoms);
	/*printf("%d %d %d\n",natoms,nchains,chainlength);*/
	if (*natoms != (nchains*chainlength))
	{
	    printf("Wrong chainlength or no_chains!\n");
	    exit(0);
	}

	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: BOX BOUNDS pp pp pp\n") != 0)
	{
		return(-3);
	}
	for(i=0; i<3; i++)
	{
		fscanf(fp, "%lg %lg\n", &box[i][0], &box[i][1]);
	}
	
	fgets(line, 1024, fp);
	/*printf("%s\n", line);*/
	if ( strcmp(line, "ITEM: ATOMS id type x y z ix iy iz\n") != 0)
	{
		return(-4);
	}
	
	return(0);
}

void printconfig(double **positions, int atoms)
{
	int atomno;
	
	for(atomno=0; atomno<atoms; atomno++)
	{
		printf("Atom %d: x=%g y=%g z=%g\n",
		       atomno+1, positions[atomno][0], positions[atomno][1], positions[atomno][2]);
	}
	
	return;
}

int readatom(double *atom, double xbox, double ybox, double zbox, int *atomtype)
{
	int atomno, molid, image[3];
	double x, y, z;


	fscanf(fp, "%d %d %lg %lg %lg %d %d %d\n",
	       &atomno, atomtype, &x, &y, &z, &image[0], &image[1], &image[2]);
	atom[0]=x+xbox*image[0];
	atom[1]=y+ybox*image[1];
	atom[2]=z+zbox*image[2];
	atom[3]=image[0];
	atom[4]=image[1];
	atom[5]=image[2];
	
	return(atomno);
}

void chaininpore(double **positions, int chainlength, int nchains, int firsttimestep, int timestep, int saveevery)
{
	int chain, monomer, atomno, inpore;
	FILE *chaintraj[nchains];

	for(chain=0; chain<nchains; chain++)
	{
		inpore=0;
		for(monomer=0; monomer<chainlength; monomer++)
		{
			atomno=chain*chainlength+monomer;
			if(positions[atomno][2] < 0)
			{
				inpore++;
			}	
		}
	
	 	char filename[20];	    
		sprintf(filename, "chain%03d.dat", chain);
		chaintraj[chain] = fopen(filename, "a");

		fprintf(chaintraj[chain], "%d %d\n", firsttimestep+timestep*saveevery, inpore);
		fclose(chaintraj[chain]);
	}   
}

double calcatomdist(double **positions){
	double delx, dely, delz, distsq, dist;
	int atom1, atom2;
	
	// just specify the atom numbers here, remembering that they should be on different chains
	atom1 = 50;
	atom2 = 250;

	delx = positions[atom1][0] - positions[atom2][0];
	dely = positions[atom1][1] - positions[atom2][1];
	delz = positions[atom1][2] - positions[atom2][2];

	distsq = delx*delx + dely*dely + delz*delz;
	dist = sqrt(distsq);

	return(dist);
}

// function to calculate whether chain ends are at bottom of pore

void calcbottom(double **positions, int chainlength, int nchains, double depth, double bdiskheight, double rho){

	// function to calculate whether chain ends are preferentially located near the pore bottom 
	
	int chain, monomer, atomno, atomid;
	int count, endcount, chaincount;
	int i, j, chains[nchains];
	double bottom, disktop;

	bottom = -1*depth;
	disktop = bottom + bdiskheight;

	endcount = 0;
	count = 0;

	for(i = 0; i < nchains; i++)
	{
		chains[i] = 0;
	}
	
	for(chain = 0; chain < nchains; chain++)
	{
		for(monomer = 0; monomer < chainlength; monomer++)
		{
			atomno = chain*chainlength + monomer;
			
			if(positions[atomno][2] < 0.0)
			{
				count++;
				
				if(chains[chain] < 1)
				{
					chains[chain] += 1;
				}
				else
				{
					continue;
				}
			}
		}	
	}

	chaincount = 0;

	for(j = 0; j < nchains; j++)
	{
		chaincount += chains[j];
	}

	for(chain = 0; chain < nchains; chain++)
	{
		for(monomer = 0; monomer < chainlength; monomer++)
		{
			atomno = chain*chainlength + monomer;

			if(positions[atomno][2] > bottom && positions[atomno][2] < disktop)
			{
				atomid = atomno + 1;
				
				if(atomid % chainlength == 0 || atomid % chainlength == 1)
				{
					endcount++;
				}
				else
				{
					continue;
				}

			}
		}
	}

	fprintf(out12, "%10.3f %d %d %d\n", rho, count, endcount, chaincount); 

}	


// function to calculate whether chains are clustered with like-chains 

void calcopening(double **positions, int chainlength, int nchains, double edgelength, double height, double diskrad, double diskheight, int firsttimestep, int timestep, int saveevery){
	
	// function to calculate which chains and chain-ends are near pore opening 
	
	int chain, monomer, atomno, atomid;
	int natoms, nends;
	int count, endcount;
	int currentstep;
	double diskvol, opening, xcenter, ycenter;
	double x, y;
	double delx, dely, rsq, r;
	double density, normdensity, enddensity, normenddensity;
	double vol, gdensity, genddensity;

	diskvol = M_PI*diskrad*diskrad*diskheight;

	// global densities of monomers and chain ends
	vol = edgelength*edgelength*height;
       	
	natoms = nchains*chainlength;
	nends = nchains*2;
	
	gdensity = natoms/vol;
	genddensity = nends/vol;	

	// z-position of the pore opening and xy coordinates of center of pore 
	opening = 0.0;
	xcenter = 0.0;
	ycenter = 0.0;

	// initialize counts and densities 
	count = 0;
	endcount = 0;
	density = 0.0;
	normdensity = 0.0;
	enddensity = 0.0;
	normenddensity = 0.0;

	for(chain = 0; chain < nchains; chain++)
	{
		for(monomer = 0; monomer < chainlength; monomer++)
		{
			atomno = chain*chainlength + monomer;

			if(positions[atomno][2] < diskheight && positions[atomno][2] > opening)
			{
				// wrap the x- and y-coordinates of the atoms 
				x = positions[atomno][0] - edgelength*positions[atomno][3];
				y = positions[atomno][1] - edgelength*positions[atomno][4];

				delx = x - xcenter;
				dely = y - ycenter;

				rsq = delx*delx + dely*dely;
				r = sqrt(rsq);

				if(r < diskrad)
				{
					count++;
					density += 1.0/diskvol;
					normdensity += (1.0/diskvol)/gdensity;

					atomid = atomno + 1;
					if(atomid % chainlength == 1)
					{
						//printf("Atom ID of chain end: %d\n", atomid); 
						endcount++;
						enddensity += 1.0/diskvol;
						normenddensity += (1.0/diskvol)/genddensity;
					}
					else if(atomid % chainlength == 0)
					{
						//printf("Atom ID of chain end: %d\n", atomid);
						endcount++;
						enddensity += 1.0/diskvol;
						normenddensity += (1.0/diskvol)/genddensity;
					}
					else
					{
						continue;
					}
				}
			}
		}
	}
	
	// print the counts and densities to the output file
	currentstep = firsttimestep + timestep*saveevery;
	fprintf(out11, "%3d %14.10f %14.10f %3d %14.10f %14.10f\n", count, density, normdensity, endcount, enddensity, normenddensity);

}

void calcbindensity(double **positions, int chainlength, int nchains, double edgelength, double binlength){

	int i, j;
	int nxbins, nybins;
	int xbin, ybin;
	int chain, monomer, atomno;
	double x, y;
	double binvol;

	// initialize the histogram bins 
	nxbins = (int) edgelength/binlength;
	nybins = (int) edgelength/binlength;

	double bindensities[nxbins][nybins], binenergies[nxbins][nybins];

	for(i = 0; i < nxbins; i++)
	{
		for(j = 0; j < nybins; j++)
		{
			bindensities[i][j] = 0.0;
			//printf("%14.10f\n",bindensities[i][j]);
			binenergies[i][j] = 0.0;
		}
	}

	binvol = binlength*binlength*binlength;
	//printf("%14.10f\n", binvol);

	for(chain = 0; chain < nchains; chain++)
	{
		for(monomer = 0; monomer < chainlength; monomer++)
		{
			atomno = chain*chainlength + monomer;

			if(positions[atomno][2] < binlength && positions[atomno][2] > 0.0)
			{
				// need to wrap the atom coordinates back into simulation box 
				x = positions[atomno][0] - edgelength*positions[atomno][3];
				y = positions[atomno][1] - edgelength*positions[atomno][4];	

				xbin = (int) floor( (-1*((-1*edgelength/2) - x))/binlength );
				ybin = (int) floor( (-1*((-1*edgelength/2) - y))/binlength );
				
				//printf("(%d, %d, %d, %14.10f, %14.10f, %14.10f, %14.10f)\n", atomno, xbin, ybin, positions[atomno][0], positions[atomno][1], x, y);

				bindensities[xbin][ybin] += 1.0/binvol;
				//printf("%14.10f\n", bindensities[xbin][ybin]);
				binenergies[xbin][ybin] += positions[atomno][6];
			}

		}
	}

	// print the densities and energies to the corresponding output files 
	for(i = 0; i < nxbins; i++)
	{
		for(j = 0; j < nybins; j++)
		{
			fprintf(out9, "%14.10f ", bindensities[i][j]);
			fprintf(out10, "%14.10f ", binenergies[i][j]); 
		}
	}

	fprintf(out9, "\n");
	fprintf(out10, "\n");
}

void calczdist(double **positions, int chainlength, int nchains, double edgelength, double height, int firsttimestep, int timestep, int saveevery, double binwidth, double bincut){
	
	int chain, monomer, atomno, natoms;
	int i, j, k, n, nbins, bin;
	double center, binvol, vol, density;
	double dist, zlo;
	char filename[20];
	
	/*FILE *distfile;
	
	sprintf(filename, "gz%010d.dat", firsttimestep+timestep*saveevery);
	distfile = fopen(filename, "w");*/

	/* initialize the histogram bins */
	nbins = (int) bincut/binwidth;
	double bincenters[nbins], bincounts[nbins], bindensity[nbins];

	for(i=0; i<nbins; i++)
	{
		bincounts[i] = 0.0;
		bindensity[i] = 0.0;
	}
	
	binvol = binwidth*edgelength*edgelength;
	natoms = nchains*chainlength;
	vol = edgelength*edgelength*height;
	density = natoms/vol;

	/* calculate and populate array with locations of bin centers */
	for(j=0; j < nbins; j++)
	{
		center = j*binwidth + binwidth/2.0;
		bincenters[j] = center;
	}	

	/* calculate the distance of the particle from the wall if within the g(z) cutoff */
	zlo = 0.0;
	for(chain=0; chain<nchains; chain++)
	{
		for(monomer=0; monomer<chainlength; monomer++)
		{
			atomno = chain*chainlength+monomer;
			dist = positions[atomno][2] - zlo;

			if(dist < bincut && dist > 0.0)
			{
				// determine which bin to add particle to
				bin = (int) floor(dist/binwidth);
				bincounts[bin] += 1.0/(binvol*density);
			        bindensity[bin] += 1.0/binvol;	
			}
		}	
	}

	// print to g(z) vs. time file initialized in main loop
	
	for(n=0; n < nbins; n++)
	{
		fprintf(out8, "%14.10f ", bindensity[n]);
	}
	fprintf(out8, "\n");
	
	// write the g(z) file
	// use this if you want to output a g(z) file periodically 
	/*fprintf(distfile, "# bincenter g(z)\n");
	for(k = 0; k < nbins; k++)
	{ 
		fprintf(distfile, "%14.10f %14.10f\n", bincenters[k], bincounts[k]);
	}
	fclose(distfile);*/

}

void calcraddist(double **positions, int chainlength, int nchains, double radius, int firsttimestep, int timestep, int saveevery){
	
	int chain, monomer, atomno;
	double c1, c2, rc;
	double delx, dely, dist, distsq, r;
	char filename[20];

	FILE *distfile;

	sprintf(filename, "raddist%09d.dat", firsttimestep+timestep*saveevery);
	distfile = fopen(filename, "w");

	/* define the center of the cylinder and the potential cutoff*/	
	c1 = 0.0;
	c2 = 0.0;
	rc = 2.5;

	for(chain=0; chain<nchains; chain++)
	{
		for(monomer=0; monomer<chainlength; monomer++)
		{
			atomno = chain*chainlength+monomer;

			delx = positions[atomno][0] - c1;
		       	dely = positions[atomno][1] - c2;
			
			distsq = delx*delx + dely*dely;
			dist = sqrt(distsq);
			r = radius - dist;

			if(r < rc && r > 0.0)
			{
				fprintf(distfile, "%10.14f\n", r);
			}	
		}
	}
	
	fclose(distfile);

}

double calcforce(double **positions, int chainlength, int nchains, double favecomp[3])
{	
	int chain, i, atomno;
	double fxcomp, fycomp, fzcomp, fxsum, fysum, fzsum;
	double fxsquared, fysquared, fzsquared, fsquared, fchain, ftotal, fave, fcomp[3];

	ftotal=0;
	fcomp[0]=0;
	fcomp[1]=0;
	fcomp[2]=0;

	for(chain=0; chain<nchains; chain++)
	{
		
		fxsum=0;
		fysum=0;
		fzsum=0;

		for(i=0; i<chainlength; i++)
		{
			atomno=i+chainlength*chain;

			
			fxcomp=positions[atomno][3];
			fycomp=positions[atomno][4];
			fzcomp=positions[atomno][5];

			fxsum+=fxcomp;
			fysum+=fycomp;
			fzsum+=fzcomp;
		}

		fcomp[0]+=fxsum;
		fcomp[1]+=fysum;
		fcomp[2]+=fzsum;

		fxsquared=fxsum*fxsum;
		fysquared=fysum*fysum;
		fzsquared=fzsum*fzsum;		

		fsquared=fxsquared+fysquared+fzsquared;

		fchain=sqrt(fsquared);
		
		ftotal+=fchain;
	}

	fave=ftotal/nchains;
	favecomp[0]=fcomp[0]/nchains;
	favecomp[1]=fcomp[1]/nchains;
	favecomp[2]=fcomp[2]/nchains;

	return(fave);
}

double calcgyration(double ***cm_array, double **positions, int chainlength, int nchains, int timestep, double coordaveRg[3])
/*Referece: polymer physics; Rubinstein p.60*/
{
	int chain, i, atomno; 
	double aveRg, x, y, z, xsum, ysum, zsum, Rg, Rgsum, coordRgsum[3];

	Rgsum=0;
	coordRgsum[0]=0;
	coordRgsum[1]=0;
	coordRgsum[2]=0;

	for(chain=0; chain<nchains; chain++)
	{
		xsum=0;
		ysum=0;
		zsum=0;
		
		for(i=0; i<chainlength; i++)
		{
			atomno=i+chainlength*chain;
			
			x=positions[atomno][0]-cm_array[timestep][chain][0];
			y=positions[atomno][1]-cm_array[timestep][chain][1];
			z=positions[atomno][2]-cm_array[timestep][chain][2];
			
			xsum+=x*x;
			ysum+=y*y;
			zsum+=z*z;
		}

		Rg=(xsum+ysum+zsum)/chainlength;
		/*printf("Chain %d's Radius of gyration: %14.10f\n", chain+1, Rg);*/

		coordRgsum[0]+=xsum/chainlength;
		coordRgsum[1]+=ysum/chainlength;
		coordRgsum[2]+=zsum/chainlength;

		Rgsum+=Rg;
		/*printf("%g %g\n",Rgsum, Rg);*/
		
	}

	coordaveRg[0]=coordRgsum[0]/nchains;
	coordaveRg[1]=coordRgsum[1]/nchains;
	coordaveRg[2]=coordRgsum[2]/nchains;

	aveRg=Rgsum/nchains;
	/*printf("Average Radius of gyration: %14.10f\n", aveRg);*/
	
	return(aveRg);
}

void calccentermass(double ***cm_array, double **positions, int chainlength, int nchains, int timestep)
{
	int chain;
	void calccm(double *cm, double **positions, int chainlength);
	
	for(chain=0; chain<nchains; chain++)
	{
		calccm(cm_array[timestep][chain], positions+chain*chainlength, chainlength);

		/*printf("chain %d's center of mass: (%14.10f, %14.10f, %14.10f)\n", 
			chain+1, cm_array[timestep][chain][0], cm_array[timestep][chain][1], cm_array[timestep][chain][2]);*/
	}

	return;
}

void calccm(double cm[3], double **positions, int chainlength)
{
	int atomno;
	double xsum, ysum, zsum;
	xsum=0;
	ysum=0;
	zsum=0;
	
	for(atomno=0; atomno<chainlength; atomno++)
	{	
		xsum+=positions[atomno][0];
		ysum+=positions[atomno][1];
		zsum+=positions[atomno][2];
	}
	cm[0]=xsum/chainlength;
	cm[1]=ysum/chainlength;
	cm[2]=zsum/chainlength;

	return;
}

double calcnumdensity(double **positions, int chainlength, int nchains, double edgelength, double height, double depth, double radius, double partcoeff, double numdensity[4], int monomers[2])
{
	int i;
	double nomon_pore, nomon_bulk;
	double porevol, bulkvol;
	double monomervol, monomerrad;
	double rad, min, offset;
	double zcut, funnelvol, fdepth;

	nomon_pore=0;
	nomon_bulk=0;

	// determine the effective pore radius by subtracting the particle diameter from 
	// the LJ-93 potential minimum and then subtract that from the value of the radius
	// set by the region command 
	min = 0.858374;
	monomerrad = 0.5;
	funnelvol = 246.983; 
	offset = min - monomerrad;
	rad = radius - offset;
	zcut = -1*radius;
	fdepth = depth - radius;

	bulkvol=edgelength*edgelength*height + funnelvol;
	porevol=M_PI*rad*rad*fdepth;

	monomervol=4*M_PI*monomerrad*monomerrad*monomerrad/3;
		
	for(i=0; i<chainlength*nchains; i++)
	{
		if (positions[i][2] < zcut)
		{
			nomon_pore++;
		}
	
		else if (positions[i][2] > zcut)
		{
			nomon_bulk++;
		}
	}

	/*printf("Number of monomers in the pore: %14.10f\n", nomon_pore);
	printf("Number of monomers in the bulk: %14.10f\n", nomon_bulk);
	printf("Volume of the pore: %10.14f\n", porevol);
	printf("Volume of the bulk: %10.14f\n", bulkvol);*/
	
	numdensity[0]=nomon_pore/porevol;
	numdensity[1]=nomon_bulk/bulkvol;

	partcoeff = numdensity[0]/numdensity[1];
	
	monomers[0]=(int)nomon_pore;

	/*printf("Number density in the bulk: %10.14f\n", numdensity[0]);
	printf("Number density in the pore: %10.14f\n", numdensitypore);*/

	return(partcoeff);
}	

double calcendlength(double **positions, int chainlength, int nchains, double aveRe, double coordaveRe[3])
/* calculate the mean square end-to-end distance of the chains*/
{
	double xend, yend, zend, Re, Resum, coordResum[3];
	int chain, n1, n2;
    

	coordResum[0]=0;
	coordResum[1]=0;
	coordResum[2]=0;
	
	Resum=0;
	
	for(chain=0; chain<nchains; chain++)
	{
		n1=chainlength*chain+chainlength-1;
		n2=chainlength*chain;	
		
		xend=positions[n1][0]-positions[n2][0];
		yend=positions[n1][1]-positions[n2][1];
		zend=positions[n1][2]-positions[n2][2];

		Re=xend*xend+yend*yend+zend*zend;
		
		/*printf("Chain %d's end-to-end length: %14.10f\n", chain+1, Re);*/

		coordResum[0]+=xend*xend;
		coordResum[1]+=yend*yend;
		coordResum[2]+=zend*zend;

		Resum+=Re;
	}

	aveRe=Resum/nchains;

	coordaveRe[0]=coordResum[0]/nchains;
	coordaveRe[1]=coordResum[1]/nchains;
	coordaveRe[2]=coordResum[2]/nchains;

	return(aveRe);
}	
double calcmonomerdistance(double **positions, int chainlength, int nchains, double aveavel, double coordaveavel[3])
{
    double xdist, ydist, zdist;
    double coordl[3], coordavel[3], distancel, suml, avel, sumcoordavel[3], sumavel;
    int chainno, atomno;

    sumcoordavel[0]=0;
    sumcoordavel[1]=0;
    sumcoordavel[2]=0;
    sumavel=0;
    
    for(chainno=0; chainno<nchains; chainno++)
    {
	coordl[0]=0;
	coordl[1]=0;
	coordl[2]=0;
	suml=0;

	for(atomno=0; atomno<(chainlength-1); atomno++)
	{
	    xdist=positions[chainno*chainlength+atomno+1][0]-positions[chainno*chainlength+atomno][0];
	    ydist=positions[chainno*chainlength+atomno+1][1]-positions[chainno*chainlength+atomno][1];
	    zdist=positions[chainno*chainlength+atomno+1][2]-positions[chainno*chainlength+atomno][2];

	    xdist=xdist*xdist;
	    ydist=ydist*ydist;
	    zdist=zdist*zdist;

	    distancel=xdist+ydist+zdist;

	    coordl[0]+=xdist;
	    coordl[1]+=ydist;
	    coordl[2]+=zdist;
	    suml+=distancel;

	}

	coordavel[0]=coordl[0]/(chainlength-1);
	coordavel[1]=coordl[1]/(chainlength-1);
	coordavel[2]=coordl[2]/(chainlength-1);
	avel=suml/(chainlength-1);

	sumcoordavel[0]+=coordavel[0];
	sumcoordavel[1]+=coordavel[1];
	sumcoordavel[2]+=coordavel[2];
	sumavel+=avel;

    }

    coordaveavel[0]=sumcoordavel[0]/nchains;
    coordaveavel[1]=sumcoordavel[1]/nchains;
    coordaveavel[2]=sumcoordavel[2]/nchains;
    aveavel=sumavel/nchains;

    return(aveavel);
}

int adsorptioncondition(double **positions, int chainlength, int nchains, int total_nomon_ad, double ad_condition, int *total_nochain_ad)
{
    int chainno, atomno, nomon_ad, nochain_ad;

    total_nomon_ad=0;
    *total_nochain_ad=0;

    for(chainno=0; chainno<nchains; chainno++)
    {
	nomon_ad=0;
	nochain_ad=0;

	for(atomno=0; atomno<chainlength; atomno++)
	{
	    if (positions[chainno*chainlength+atomno][2] <= ad_condition)
	    {
		nomon_ad++;
	    }
	}

	if (nomon_ad >= 1)
	{
	    nochain_ad++;
	}
	total_nomon_ad+=nomon_ad;
	*total_nochain_ad+=nochain_ad;
    }

    /*printf("%d\n", *total_nochain_ad);*/
    return(total_nomon_ad);
}

void looptraintail(double **positions, int chainlength, int nchains, double ad_condition, int *no_train_mon, int *no_tail_mon, int *no_loop_mon, int *no_nonadsorb_mon, double *avetrain_mon, double *avetail_mon, double *aveloop_mon, double *avenonadsorb_mon, int *no_train, int *no_tail, int *no_loop, int *no_nonadsorb)
{

    int chainno, atomno, index;
    int monomer[chainlength*nchains];
    int train[chainlength*nchains], tail[chainlength*nchains], loop[chainlength*nchains], nonadsorb[chainlength*nchains];
    int a, b, c, d;
        
    for(chainno=0; chainno<nchains; chainno++)
    {
	for(atomno=0; atomno<chainlength; atomno++)
	{
	    index=chainno*chainlength+atomno;

	    if (positions[index][2] > ad_condition)
	    {
		monomer[index]=0;

		if (monomer[index-1]==1)
		{
		    if (index == chainno*chainlength)
			goto end;
	
		    monomer[index]=2;
		}
		else if (monomer[index-1]==2)
		{
		    if (index == chainno*chainlength)
			goto end;	

		    monomer[index]=2;
		}
	    }

	    else if (positions[index][2] <= ad_condition)
	    {
		monomer[index]=1;
		
		if (monomer[index-1]==0)
		{
		    while (monomer[index-1]==0)
		    {
			if(index == chainno*chainlength)
			    goto end;
				
			monomer[index-1]=2;
			index--;
		    }
		}
		else if (monomer[index-1]==2)
		{
		    while (monomer[index-1]==2)
		    {
			if(index == chainno*chainlength)
			    goto end;
					
			monomer[index-1]=3;
			index--;
		    }
		}
	    }
	end: ;
	}
    }  

    *no_train=*no_tail=*no_loop=*no_nonadsorb=0;
    
    for(chainno=0; chainno<nchains; chainno++)
    {
    	for(atomno=0; atomno<chainlength; atomno++)
	{ 
	    index=chainno*chainlength+atomno;
	    
	    switch(monomer[index])
	    {
		case(0):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit0;
			
			nonadsorb[*no_nonadsorb]+=1;
		    }
		    else
		    {
		    exit0:
		
			*no_nonadsorb+=1;
			nonadsorb[*no_nonadsorb]=1;
		    }
		    break;

		case(1):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit1;
			
			train[*no_train]+=1;
		    }
	 
		    else 
		    {
		    exit1:
			
			*no_train+=1;
			train[*no_train]=1;
		    }
		    break;

		case(2):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit2;
			
			tail[*no_tail]+=1;
		    }
		    else
		    {
		    exit2:
			
			*no_tail+=1;
			tail[*no_tail]=1;
		    }
		    break;

		case(3):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit3;
			
			loop[*no_loop]+=1;
		    }
		    else
		    {
		    exit3:
			
			*no_loop+=1;
			loop[*no_loop]=1;
		    }
		    break;
	    }
	}
    }

    /*printf("%d %d %d %d\n", no_nonadsorb, no_train, no_tail, no_loop);*/
   
    *no_nonadsorb_mon=*no_train_mon=*no_tail_mon=*no_loop_mon=0;
    *avenonadsorb_mon=*avetrain_mon=*avetail_mon=*aveloop_mon=0;
   
    for(a=0; a<*no_nonadsorb; a++)
    {
	*no_nonadsorb_mon+=nonadsorb[a+1];
    }
    if(*no_nonadsorb != 0)
	*avenonadsorb_mon=(float) *no_nonadsorb_mon/ *no_nonadsorb;
 
    for(b=0; b<*no_train; b++)
    {
	*no_train_mon+=train[b+1];
    }
    if(*no_train != 0)
    *avetrain_mon=(float) *no_train_mon/ *no_train;
    
    for(c=0; c<*no_tail; c++)
    {
	*no_tail_mon+=tail[c+1];
    }
    if(*no_tail != 0)
    *avetail_mon=(float) *no_tail_mon/ *no_tail;

    for(d=0; d<*no_loop; d++)
    {
	*no_loop_mon+=loop[d+1];
    }
    if(*no_loop != 0)
    *aveloop_mon=(float) *no_loop_mon/ *no_loop;

    /*printf("%d %d %d %d\n", no_nonadsorb_mon, no_train_mon, no_tail_mon, no_loop_mon);*/

     /* printf("%10.4f %10.4f %10.4f %10.4f\n", avenonadsorb_mon, avetrain_mon, avetail_mon, aveloop_mon);*/

    return;
}
