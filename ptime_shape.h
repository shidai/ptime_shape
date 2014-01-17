#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
#include <fftw3.h>
#include "fitsio.h"
//#include "ptime.h"

#define NP 2048
#define PI 3.14159265359

long int stt_imjd ( char *name );
double read_psrfreq ( char *name );

int get_nchan ( char *name );
int get_npol ( char *name );
int get_nphase ( char *name );
int get_subint ( char *name );

int read_scl ( char *name, int subint, double *scl, int nchan, int npol);
int read_offs ( char *name, int subint, double *offs, int nchan, int npol);
int read_value ( char *name, int subint, double *value, int nphase, int nchan, int npol);

int read_prof ( char *name, int subint, double *profile, int nphase, int npol, int nchan);

int check_std ( char *name, int subint, int mode, int nchn, int nphase, int npol);

int read_std_pt ( char *name, double *profile, int nphase, int nchn, int npol);

int read_std ( char *name, int subint, double *profile, int nphase, int mode, int nchn, int npol);

int find_peak (int n, double *s, int *position);
double find_peak_value (int n, double *s);

int on_pulse (int nphase, int peak_position, double *in, double *out, double frac);
int off_pulse (int nphase, int peak_position, double *in, double *out, double frac_off, double frac_on);

int remove_baseline (double *in, double frac_off, double frac_on, int n, double *out);

int error (double *p_off, double *s_on, double scale, double ro, double snr, int num_on, int num_off, double *err);
int shape_para (double *s, double *p, int nphase, double frac_on, double frac_off, FILE *fp, double psrfreq, long int mjd, int nchn, int npol, int nsub, char *fname);

int real_obs (char *fname, char *tname, char *oname, int mode, FILE *fp, double frac_on, double frac_off);

//int simulate (int n, double SNR, double *s, double *p);
//int sim_obs (char *tname, char *oname, int mode, double SNR, FILE *fp, double psrfreq);

/////////////////////////////////////////////////////////////////
//
double A7 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

double A9 (double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int dft_profiles (int N, double *in, fftw_complex *out);

double zbrent(double (*func)(double phase, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn), double x1, double x2, double tol, double a_s[][NP], double a_p[][NP], double p_s[][NP], double p_p[][NP], int num, int nchn);

int preA7 (int *k, double amp_s[][NP], double amp_p[][NP], double phi_s[][NP], double phi_p[][NP], double *s, double *p, int nphase, int nchn, double *real_p, double *ima_p);

int inverse_dft (double *real_p, double *ima_p, int ncount, double *p_new);

int align (int N, double phase, double b, double a, double *real_p, double *real_p_align, double *ima_p, double *ima_p_align);

int get_toa (double *s, double *p, double *p_new, double *scale, double psrfreq, int nphase, long int mjd, int nchn, int npol, int nsub, char *fname);
