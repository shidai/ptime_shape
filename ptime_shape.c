//"Usage: ptime_shape -f fname -std tname (-pt tname) -o oname \n"
//"Calculate the shape parameter\n"
//fname: data file; tname: templates; oname: output .tim; -std: standard template format; -pt: ptime template;\n"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ptime_shape.h"

int main (int argc, char *argv[])
{
	//int h,i,j,k;

	//////////////////////////////////////////////////////
	char fname[128];   // name of data file
	char tname[128];   // name of template
	char oname[128];   // name of output .tim
	int mode; // to distinguish different type of templates
	int smode = 0; // if smode==1, simulate profiles

	double frac_on, frac_off;

	int i;
	int index, n;
	for (i=0;i<argc;i++)
    {
		if (strcmp(argv[i],"-f") == 0)
		{
            index = i + 1;
			n = 0;
			while ( (index + n) < argc && strcmp(argv[index+n],"-std") != 0 && strcmp(argv[index+n],"-pt") != 0 && strcmp(argv[index+n],"-o") != 0 && strcmp(argv[index+n],"-sim") != 0 && strcmp(argv[index+n],"-frac_on") != 0 && strcmp(argv[index+n],"-frac_off") != 0)
			{
				n++;
		    }
			//strcpy(fname,argv[++i]);
		}
		else if (strcmp(argv[i],"-std")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 0; // standard template format
			printf ("standard template format\n");
			//sscanf(argv[++i],"%d",&nbin);
		}
		else if (strcmp(argv[i],"-pt")==0)
		{
			strcpy(tname,argv[++i]);
			mode = 1; // ptime template
			printf ("ptime template format\n");
		}
		else if (strcmp(argv[i],"-o")==0)
		{
			strcpy(oname,argv[++i]);
		}
		else if (strcmp(argv[i],"-sim")==0)
		{
			smode = 1;
		}
		else if (strcmp(argv[i],"-frac_on")==0)
		{
			frac_on = atof(argv[++i]);
		}
		else if (strcmp(argv[i],"-frac_off")==0)
		{
			frac_off = atof(argv[++i]);
		}
    }
	//printf ("%d\n", smode);

	// start to deal with different data file
	//
	// open file to write toa 
	FILE *fp;
	if ((fp = fopen(oname, "w+")) == NULL)
	{
        fprintf (stdout, "Can't open file\n");
		exit(1);
	}
    //fprintf (fp, "S0    S    err\n");
	/////////////////////////////////////////////////////////
	
	if (smode == 0)
	{
		int k;
		for (k = index; k < index + n; k++)
		{
			// get the data file name
			strcpy(fname,argv[k]);
			real_obs(fname, tname, oname, mode, fp, frac_on, frac_off);
		}
	}
	else if (smode == 1)
	{
		double SNR;
		for (SNR = 100.0; SNR < 1000.0; SNR += 100.0)
		{
			printf ("SNR: %lf\n", SNR);
			sim_obs (tname, oname, mode, SNR, fp, 1000);
		}
	}

    if (fclose (fp) != 0)
		fprintf (stderr, "Error closing\n");

	return 0;
}
