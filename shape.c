// calculate the stability of profile
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include "ptime_shape.h"

int find_peak (int n, double *s, int *position)
{
	int i;
	double temp[n];
	double peak;

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}
	peak=temp[n-1];

	for (i=0;i<n;i++)
	{
		if (fabs(peak-s[i])<1.0e-3)
		{
			(*position)=i;
		}
	}

	return 0;
}

double find_peak_value (int n, double *s)
{
	int i;
	double temp[n];

	for (i=0;i<n;i++)
	{
		temp[i]=s[i];
	}

	double a,b,c;
	for (i=0;i<n-1;i++)
	{
		a=temp[i];
		b=temp[i+1];
		c=(a>=b ? a : b);

		temp[i+1]=c;
	}

	return temp[n-1];
}

int on_pulse (int nphase, int peak_position, double *in, double *out, double frac)
// define the on_pulse range, 50% of the phase
{
	int n = nphase;
	int num = (int)(n*frac);
	int i;
	for (i = 0; i < num; i++)
	{
		if ((peak_position-num/2+i) < 0)
		{
			out[i] = in[(n-1)+(peak_position-num/2+i)];
		}
		else if ((peak_position-num/2+i) > n-1)
		{
			out[i] = in[(peak_position-num/2+i)-(n-1)];
		}
		else
		{
			out[i] = in[peak_position-num/2+i];
		}
	}

	return 0;
}

int off_pulse (int nphase, double *in, double *out, double frac)
// define the on_pulse range, 10% of the phase
{
	int n = nphase;
	int num = (int)(n*frac);
	int i,j;
	double small;
	double ave;
	int index = 0;

	for (i = 0; i < n; i++)
	{
		if (i == 0)
		{
			small = 0.0;
			for(j = 0; j < num; j++)
			{
				small += in[j];
			}
			small = small/num;
		}
			
		ave = 0.0;
		for(j = 0; j < num; j++)
		{
			if ((i+j) > n-1)
			{
				ave += in[(i+j)-(n-1)];
			}
			else 
			{
				ave += in[i+j];
			}
		}
		ave = ave/num;

		small = (ave <= small ? ave : small);
		index = (ave <= small ? i : index);
		//printf ("%d %lf %lf\n", index, small, ave);
	}

	for (i = 0; i < num; i++)
	{
		if ((index+i) > n-1)
		{
			out[i] = in[(index+i)-(n-1)];
		}
		else 
		{
			out[i] = in[index+i];
		}
	}

	return 0;
}

int remove_baseline (double *in, double frac_off, int n, double *out)
{
	// define the off_pulse range of std, frac_off is the fraction of the phase
	int num_off = (int)(n*frac_off);
	double off_0[num_off];

	off_pulse(n, in, off_0, frac_off);

	// normalize std
	int i;
	double baseline = 0.0;
    for (i = 0; i < num_off; i++)
    {
        baseline += off_0[i];
        //average_s += s_off[i];
    }
	baseline = baseline/num_off;

    printf ("the baseline of std is: %lf \n", baseline);
    //printf ("average is: %lf %lf\n", average, average_s);

	for (i = 0; i < n; i++)
	{
		out[i] = (in[i]-baseline);
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	return 0;
}

int error (double *p_off, double *s_on, int num_on, int num_off, double *err)
{
	int i;
    double s_bar=0.0;
    
    for (i = 0; i < num_on; i++)
    {
        s_bar += s_on[i];
    }

    s_bar=s_bar/num_on;

    int m=0;
	double average = 0.0;
	//double average_s = 0.0;
    for (i = 0; i < num_off; i++)
    {
        average += p_off[i];
        //average_s += s_off[i];
		m++;
    }
	average = average/m;
	//average_s = average_s/m;

    double sigma=0.0;
    for (i = 0; i < num_off; i++)
    {
        //rms += (p_off[i])*(p_off[i]);
        sigma += (p_off[i]-average)*(p_off[i]-average);
    }

    sigma=sqrt(sigma/m);

	double ro21;
	ro21 = 0.0;
    for (i = 0; i < num_on; i++)
    {
        ro21 += (s_on[i] - s_bar)*(s_on[i] - s_bar);
    }

	double chi;
	chi = sigma*sigma/ro21;
	
	(*err) = chi*sqrt((num_on-1.0)/(2.0*num_on*num_on));

	return 0;
}

int shape_para (double *s, double *p, int nphase, double frac_on, double frac_off, FILE *fp, double psrfreq, long int mjd, int nchn, int npol)
{
	int n = nphase;
	
	// remove the baseline
	double s_nobase[n], p_nobase[n];
	remove_baseline (s, frac_off, n, s_nobase);
	remove_baseline (p, frac_off, n, p_nobase);

	// find the peak and peak position
	int s_peak_position, p_peak_position;
	double s_peak, p_peak;
	
	s_peak = find_peak_value (n, s_nobase);
	p_peak = find_peak_value (n, p_nobase);

	find_peak (n, s_nobase, &s_peak_position);
	find_peak (n, p_nobase, &p_peak_position);

	// normalize std
	int i;
	double s_norm[n];

	for (i = 0; i < n; i++)
	{
		s_norm[i] = s_nobase[i]/s_peak;
		//printf ("%d %lf %lf\n", i, s[i], s_nobase[i]);
		//s_norm[i] = s_nobase[i]/s_peak;
		//s_norm[i] = (s[i]-baseline)/(s_peak-baseline);
	}
	
	// define the off_pulse range, frac_off is the fraction of the phase
	int num_off = (int)(n*frac_off);
	double s_off[num_off], p_off[num_off];

	off_pulse(n, s_norm, s_off, frac_off);
	off_pulse(n, p_nobase, p_off, frac_off);

	// define the on_pulse range, frac_on is the fraction of the phase
	int num_on = (int)(n*frac_on);
	double s_on[num_on], p_on_0[num_on], p_on[num_on];

	on_pulse(n, s_peak_position, s_norm, s_on, frac_on);
	on_pulse(n, p_peak_position, p_nobase, p_on_0, frac_on);

	get_toa (s_on, p_on_0, p_on, psrfreq, num_on);
	//printf ("///////////////////\n");

	/*
	for (i = 0; i < num_on; i++)
	{
		printf ("%d %lf %lf %lf\n", i, s_on[i], p_on[i], p_on_0[i]);
	}
	for (i = 0; i < (int)(n/2); i++)
	{
		printf ("%d %lf\n", i, p_on[i]);
	}

    /////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	// test things
    double s_bar_t=0.0;
    double p_bar_t=0.0;
    
    //for (i = 0; i < num_on; i++)
    for (i = 0; i < num_off; i++)
    {
        //s_bar_t += s_on[i];
        //p_bar_t += p_on[i];
        s_bar_t += s_off[i];
        p_bar_t += p_off[i];
    }

    //s_bar_t = s_bar_t/num_on;
    //p_bar_t = p_bar_t/num_on;
    s_bar_t = s_bar_t/num_off;
    p_bar_t = p_bar_t/num_off;
    //fprintf (fp, "%lf %lf\n", s_bar_t, p_bar_t);

    double ss = 0.0;
    double pp = 0.0;

    for (i = 0; i < num_on; i++)
    {
        ss += (s_on[i] - s_bar_t)*(s_on[i] - s_bar_t);
        pp += (p_on[i] - p_bar_t)*(p_on[i] - p_bar_t);
		//fprintf (fp, "%lf %lf\n", s_on[i] - s_bar_t, p_on[i] - p_bar_t);
    }

	double test1, test2;
    test1 = sqrt(num_on/(2.0*ss));
    test2 = 2000*sqrt(num_on/(2.0*pp));

    fprintf (fp, "%lf %lf %lf %lf\n", test1, test2, s_bar_t, p_bar_t);
	*/
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////
    double s_bar=0.0;
    double p_bar=0.0;
    
    //for (i = 0; i < num_on; i++)
    for (i = 0; i < num_off; i++)
    {
        //printf ("s is: %f\n", s[i]);
        //printf ("p is: %f\n", p[i]);
        s_bar += s_off[i];
        p_bar += p_off[i];
        //printf ("s_bar is: %f\n", s_bar);
        //printf ("p_bar is: %f\n", p_bar);
    }

    s_bar=s_bar/num_off;
    p_bar=p_bar/num_off;
    printf ("s_bar is: %f\n", s_bar);
    printf ("p_bar is: %f\n", p_bar);

    //////////////////////////////////////////////////////////////////////////////
    double ro1=0.0;
    double ro21=0.0;
    double ro22=0.0;
    double ro;

    for (i = 0; i < num_on; i++)
    {
		//printf ("%d %lf %lf \n", i, s_on[i]-s_bar, p_on[i]-p_bar);
		ro1 += (s_on[i] - s_bar)*(p_on[i] - p_bar);
        ro21 += (s_on[i] - s_bar)*(s_on[i] - s_bar);
        ro22 += (p_on[i] - p_bar)*(p_on[i] - p_bar);
    }

    ro=ro1/(sqrt(ro21)*sqrt(ro22));
    //ro_s=1.0/sqrt(ro21);
    printf ("ro is: %lf \n", ro);
    //printf ("ro is: %lf %lf\n", ro, ro_s);
    
    ///////////////////////////////////////////////////////////////////////////////////

    double rms=0.0;
    //double rms_s=0.0;
    for (i = 0; i < num_off; i++)
    {
        //rms += (p_off[i])*(p_off[i]);
        rms += (p_off[i]-p_bar)*(p_off[i]-p_bar);
        //rms_s += (s_off[i]-average_s)*(s_off[i]-average_s);
    }

    rms=sqrt(rms/num_off);
    //rms_s=sqrt(rms_s/m);
    printf ("rms is: %lf \n", rms);
    //printf ("rms is: %lf %lf\n", rms, rms_s);

    ////////////////////////////////////////////////////////////////////////////////////

    double shape_para;
    double shape_para_std;
    //double SNR = (p_peak)/rms;
    double SNR = (p_peak-p_bar)/rms;
    //double SNR_s = (s_peak-average_s)/rms_s;
    printf ("SNR is: %lf \n", SNR);
    //printf ("SNR is: %lf %lf\n", SNR, SNR_s);

    shape_para = SNR*sqrt(1.0-ro);
    //shape_para_std = SNR_s*sqrt(1-ro_s);
    shape_para_std = sqrt(num_on/(2.0*ro21));
    double shape_para_test = p_peak*sqrt(num_on/(2.0*ro22));

	double err;
	error (p_off, s_on, num_on, num_off, &err);
    //fprintf (fp, "The shape parameter is: %f\n", shape_para);
    //fprintf (fp, "The theoretical shape parameter is: %f\n", shape_para_std);
    //fprintf (fp, "%lf %lf %.10lf\n", shape_para_std, shape_para, err);
    fprintf (fp, "%ld %d %d %.3lf %.3lf %.3lf %.10lf\n", mjd, nchn, npol, shape_para_std, shape_para, shape_para_test, err);

    return 0;
}

int real_obs (char *fname, char *tname, char *oname, int mode, FILE *fp, double frac_on, double frac_off)
{
	// name of different extension of data files
	char name_data[50]; 
	char name_psrparam[50]; 

	char data[] = "[SUBINT]";
	char psrparam[] = "[PSRPARAM]";

	// read name of std
	char std[50];
	strcpy(std,tname);
	if ( mode == 0)
	{
		strcat(std, data);
	}

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	int h, i, j, p;
	{
		// name of different extension
		strcpy(name_data,fname);
		strcpy(name_psrparam,fname);

		strcat(name_data, data);
		strcat(name_psrparam, psrparam);

		////////////////////////////////////////////////////
		
		double psrfreq;
		psrfreq = read_psrfreq(name_psrparam);
		printf ("PSR frequency: %.15lf\n", psrfreq);

		////////////////////////////////////////////////////

		long int imjd;
		imjd = stt_imjd(fname);
	
		int nphase;
		int nchn;
		int nsub;
		int npol;
	
		nchn = get_nchan(name_data);	
		npol = get_npol(name_data);	
		nsub = get_subint(name_data);	
		nphase = get_nphase(name_data);	

		//printf ("%d\n", nchn);
		////////////////////////////////////////////////

		// read a std
		double s_multi[nphase*nchn*npol];
		double s_temp[nphase];

		//readfile(argv[1],&n,tt,s);
		//read_prof(std,1,s_multi,nphase);
	
		// check the channel and phase number of template
		check_std(std,1,mode,nchn,nphase,npol);

		read_std(std,1,s_multi,nphase,mode,nchn,npol);

		////////////////////////////////////////////////////////////////////////////////

		double p_multi[nchn*npol*nphase];
		double p_temp[nphase];

		// start to calculate shape paremeter for different subint
		for (h = 1; h <= nsub; h++)
		{
			// simulate data

			//SNR = 500.0 + 200.0*i;
			//simulate(n,SNR,s,p_temp);

			// read profiles from data file
			read_prof(name_data,h,p_multi,nphase,npol,nchn);
			//readfile(argv[2],&n,tt,p_multi);

			// start to calculate shape parameter for different channels
			for (i = 0; i < nchn; i++)
			{
				for (p = 0; p < npol; p++)
				{
					for (j = 0; j < nphase; j++)
					{
						//printf ("%lf %lf\n", p_multi[j], s[j]);
						//s_multi[i*nphase + j] = s[j];
						p_temp[j] = p_multi[i*npol*nphase + p*nphase + j];
						s_temp[j] = s_multi[i*npol*nphase + p*nphase + j];
						//s_temp[j] = s_multi[i*nphase + j];
						//fprintf (fp, "%d %d %lf\n", i, j, p_temp[j]);
					}
					//get_toa (s_temp, p_temp, p_new, psrfreq, nphase);
					shape_para(s_temp, p_temp, nphase, frac_on, frac_off, fp, psrfreq, imjd, i, p);
				}
			}
		}
	}

	return 0;
}

/*
int simulate (int n, double SNR, double *s, double *p)
{
	// simulate a profile with white noise
	///////////////////////////////////////////////////////////////////////
	// initialize gsl 
	
	int i;
	const gsl_rng_type * T;
	gsl_rng * r;

	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
  
	////////////////////////////////////////////////////////////////////////
	//  determine the amplitude of white noise according to SNR
	
	double scale;   // the scale multiply to white noise to get certain SNR
	double amp_noise, noise[n];
	double peak_s;

	amp_noise = 0.0;
	for (i=0;i<n;i++)
	{
		noise[i]=gsl_ran_gaussian(r,1.0);
		//printf ("%d %lf\n", i, noise[i]);
		amp_noise+=noise[i]*noise[i];
	}
	
	amp_noise=sqrt(amp_noise/n);
    //printf ("%g\n", amp_noise);

	peak_s=find_peak_value(n,s);   // find the peak flux of the std
    //printf ("peak of std: %g\n", peak_s);

	scale=peak_s/(SNR*amp_noise);
    //printf ("%g\n", scale);
	
	//////////////////////////////////////////////////////////////////////////
	//  add noise to std ==> p

	for (i=0;i<n;i++)
	{
		p[i]=(s[i]+scale*noise[i]);
		//printf ("%d %lf\n", i, p[i]);
	}

	gsl_rng_free (r);
  
	return 0;
}

int sim_obs (char *tname, char *oname, int mode, double SNR, FILE *fp, double psrfreq)
{
	char data[] = "[SUBINT]";

	// read name of std
	char std[50];
	int nphase;
	int nchn;
	int nsub;
	int npol;

	strcpy(std,tname);

	if ( mode == 0)
	{
		strcat(std, data);

		nchn = get_nchan(std);	
		npol = get_npol(std);	
		nsub = get_subint(std);	
		nphase = get_nphase(std);	
	}
	else 
	{
		nchn = 1;	
		npol = 1;	
		nsub = 1;	
		nphase = 1024;	
	}

	/////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////
	int h, i, j, p;
	{
		////////////////////////////////////////////////////

		// read a std
		double s_multi[nphase*nchn*npol];
		double s_temp[nphase];

		//readfile(argv[1],&n,tt,s);
		//read_prof(std,1,s_multi,nphase);
	
		// check the channel and phase number of template
		//check_std(std,1,mode,nchn,nphase);

		read_std(std,1,s_multi,nphase,mode,nchn,npol);

		////////////////////////////////////////////////////////////////////////////////

		//double p_multi[nchn*npol*nphase];
		double p_temp[nphase];

		// start to calculate shape paremeter for different subint
		for (h = 1; h <= nsub; h++)
		{
			// simulate data

			//SNR = 500.0 + 200.0*i;
			//simulate(n,SNR,s,p_temp);

			// read profiles from data file
			//read_prof(name_data,h,p_multi,nphase,npol,nchn);
			//readfile(argv[2],&n,tt,p_multi);

			// start to calculate shape parameter for different channels
			for (i = 0; i < nchn; i++)
			{
				for (p = 0; p < npol; p++)
				{
					for (j = 0; j < nphase; j++)
					{
						//printf ("%lf %lf\n", p_multi[j], s[j]);
						//s_multi[i*nphase + j] = s[j];
						//s_temp[j] = s_multi[i*nphase + j];
						s_temp[j] = s_multi[i*npol*nphase + p*nphase + j];
						//fprintf (fp, "%d %d %lf\n", i, j, p_temp[j]);
					}
					simulate(nphase, SNR, s_temp, p_temp);
					shape_para(s_temp, p_temp, nphase, 0.8, 0.2, fp, psrfreq);
					//printf ("%d %d\n", nphase, nchn);
				}
			}
		}
	}

	return 0;
}
*/
