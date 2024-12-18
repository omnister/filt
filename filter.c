/*
 *  Filter.c:  Classical IIR filters
 *
 *  Written by Rick Walker, HP Labs
 *  Modifications by Scott Willingham, HP Labs
 *
 *  Mods by s.d.willingham, Thu Aug 12 11:32:21 PDT 1993
 *    - Added filter initial condition option (-i)
 *    - Revised misc. code and notations ( z=exp(jwT) !)
 *    - Revised magnitude scaling factors
 *    - Added several filter approximations: maxflat, modcheb, and invcheb
 *    - changed "dump" format to be consistent with z=exp(jwT) notation.
 *    - fixed bugs in highpass, bandpass, and band-reject transformations
 *    - recoded most of initquad() and initpole()
 *    - outputs impulse/step response in autoplot format (-I or -S)
 *    - magnitude+delay or loss+delay plots (-X, -L, -P)
 *    - fixed bug in Ztransform(), where 1st-order expressions were
 *      transformed incorrectly.  (They were formerly transformed into
 *      a 2nd-order function with an extra cancelling pole-zero at z = -1.)
 *    - fixed getquad() to recognize proper order of filter when odd
 *    - Unglobalized some variables (8/19/93)
 *    - Added switch (-Fs) to normalize all frequencies by Fs (8/19/93)
 *    - consolidated filter_type & type into one variable (8/20/93)
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "bessel.h"		/* coefficients for bessel approximation */

double T = 1.0;			/* normalized sampling time period */
double fnorm = 0.0;		/* normalization frequency [set in main()] */

#define FSAMPLE  (1.0/T)
#define FNYQUIST (0.5/T)

/* macros for filter_type manipulation */
#define UNKNOWN	    0
#define LP 	    (1<<0)
#define HP 	    (1<<1)
#define BP 	    (1<<2)
#define BR 	    (1<<3)
#define FT_BAND     (LP|HP|BP|BR)
#define BUTTERWORTH (1<<4)
#define MAXFLAT	    (1<<5)
#define CHEBYSHEV   (1<<6)
#define MOD_CHEB    (1<<7)
#define INV_CHEB    (1<<8)
#define BESSEL      (1<<9)
#define FT_APPROX   (BUTTERWORTH|MAXFLAT|CHEBYSHEV|MOD_CHEB|INV_CHEB|BESSEL)

#define PLOT_MAGN  (1<<0)
#define PLOT_LOSS  (1<<1)
#define PLOT_DELAY (1<<2)
#define PLOT_PBAND (1<<3)

#define MAXBUF 1024

struct biquad {
    double a1;
    double a2;
    double b0;
    double b1;
    double b2;
    double v1;
    double v2;
    struct biquad *next;    /* next entry in chain */
}; 

struct biquad *quadlist;
struct biquad *Ztransform(); 
struct biquad *initpole(); 
struct biquad *initquad(); 
struct biquad *qalloc();
struct biquad *getquad();
double fevalquad();
double devalquad();
double warp();
double prewarp();
double filter();
double initfilter();
double asinh();
double ChebPoly();

char *progname;
int scan_OK;
char instring[MAXBUF];


main(argc, argv)    /* nth order time-domain filter */
int argc;
char *argv[];
{
    char *progname3;
    extern int  optind;  /* argv index of next option */
    extern int  opterr;
    extern char *optarg;
    int c;
    int errflag;
    FILE *fp = NULL;

  /********** Defaults **********/
    int filter_type = BUTTERWORTH;
    int order = 8;	    	    /* default filter order */
    double cf = 0.5;		    /* default center frequency */
    double bw = 0.1;	    	    /* default bw frequency */
    double ripple = 3.0103;	    /* default db ripple */
    double sb_ratio = 2.0;	    /* default fsb/fpb ratio for invcheb */
    double dc = 0.0;		    /* default dc inital state */
    int xfer_plot=0;
    int Isamples=0;		    /* plot impulse response if != 0 */
    int Ssamples=0;		    /* plot step response if != 0 */
    int dump_flag=0;
    int quiet_flag=0;
  /********** End Defaults **********/

    double eps;		/* ripple coefficient */
    double mf_radius;	/* radius of maxflat pole locus */

    int n, nquads;
    double theta, Omegax_sq, lambda, g, u,v, x,y;
    double time, F, out, outr, outi, delay;

    quadlist = (struct biquad *) NULL;	/* initialize head of list */

    errflag = 0;
    opterr = 0;     /* disables getopt's error msg's */

    progname = argv[0];

    progname3 = progname;
    if (strlen(progname3) > 3)
	progname3 += strlen(progname3) - 3;

    if (!strcmp(progname3,"lpf")) {
	filter_type = (filter_type & ~FT_BAND) | LP;
    } else if (!strcmp(progname3,"hpf")) {
	filter_type = (filter_type & ~FT_BAND) | HP;
    } else if (!strcmp(progname3,"bpf")) {
	filter_type = (filter_type & ~FT_BAND) | BP;
    } else if (!strcmp(progname3,"brf")) {
	filter_type = (filter_type & ~FT_BAND) | BR;
    } else {
	filter_type = (filter_type & ~FT_BAND) | LP; /* default */
    }

    /* parse args */
    while ((c = getopt(argc, argv, "t:n:c:b:r:s:F:i:xXLPI:S:df:q")) != EOF)
        switch (c) {
        case 't':	    /* get filter approximation type */
 	    scan_OK = sscanf(optarg,"%s", instring);
            if (scan_OK != 1) {
		errflag++;
	    } else if (!strncmp(instring,"butterworth",2)) {
		filter_type = (filter_type & ~FT_APPROX) | BUTTERWORTH;
	    } else if (!strncmp(instring,"maxflat",2)) {
		filter_type = (filter_type & ~FT_APPROX) | MAXFLAT;
	    } else if (!strncmp(instring,"bessel",2)) {
		filter_type = (filter_type & ~FT_APPROX) | BESSEL;
	    } else if (!strncmp(instring,"chebyshev",2)) {
		filter_type = (filter_type & ~FT_APPROX) | CHEBYSHEV;
	    } else if (!strncmp(instring,"modcheb",2)) {
		filter_type = (filter_type & ~FT_APPROX) | MOD_CHEB;
	    } else if (!strncmp(instring,"invcheb",2)) {
		filter_type = (filter_type & ~FT_APPROX) | INV_CHEB;
	    } else {
        	fprintf(stderr, 
		    "%s: error: filter type must be either butterworth, maxflat,\nchebyshev, modcheb, invcheb, or bessel\n\n", progname);
		errflag++;
	    }
            break;
        case 'n':	    /* filter order */
 	    scan_OK = sscanf(optarg,"%d", &order);
            if (scan_OK != 1) errflag++;
    	    if (order<=0) {
        	fprintf(stderr, 
		    "%s: error: order must be greater than 0\n", progname);
	    }
            break;
        case 'c':	    /* get center frequency for BP, BR */
 	    scan_OK = sscanf(optarg,"%lf", &cf);
            if (scan_OK != 1) errflag++;
            break;
        case 'b':	    /* filter bw frequency */
 	    scan_OK = sscanf(optarg,"%lf", &bw);
            if (scan_OK != 1) errflag++;
            break;
        case 'r':	    /* get ripple for chebyshev, maxflat */
 	    scan_OK = sscanf(optarg,"%lf", &ripple);
            if (scan_OK != 1) errflag++;
            break;
        case 's':	    /* sband/pband ratio for invcheb */
 	    scan_OK = sscanf(optarg,"%lf", &sb_ratio);
            if (scan_OK != 1) errflag++;
            break;
        case 'F':	    /* normalize frequencies to Fs */
	    if (optarg) {
		switch (*optarg) {
		    case 'n': case 'N':  fnorm = FNYQUIST; break;
		    case 's': case 'S':  fnorm = FSAMPLE; break;
		    default:   errflag++;
		}
	    } else
		errflag++;
            break;
        case 'i':	    /* Initialize "dc state" */
 	    scan_OK = sscanf(optarg,"%lf", &dc);
            if (scan_OK != 1) errflag++;
            break;
        case 'x':	    /* don't filter, just print out xfer function */
            xfer_plot = PLOT_MAGN;
            break;
        case 'X':	    /* xfer function: magnitude + delay */
            xfer_plot = PLOT_MAGN | PLOT_DELAY;
            break;
        case 'L':	    /* loss + delay */
            xfer_plot = PLOT_LOSS | PLOT_DELAY;
            break;
        case 'P':	    /* loss + delay, passband only */
            xfer_plot = PLOT_LOSS | PLOT_DELAY | PLOT_PBAND;
            break;
        case 'I':	    /* print out impulse response */
 	    scan_OK = sscanf(optarg,"%d", &Isamples);
            if (scan_OK != 1) errflag++;
            break;
        case 'S':	    /* print out step response */
 	    scan_OK = sscanf(optarg,"%d", &Ssamples);
            if (scan_OK != 1) errflag++;
            break;
        case 'd':	    /* don't filter, just print out biquad values */
             dump_flag++;
	     break;
        case 'f':	    /* get input file name */
 	    scan_OK = sscanf(optarg,"%s", instring);
            if (scan_OK != 1) {
		errflag++;
	    }
	    if (( fp = fopen(instring, "r")) == NULL) {
        	fprintf(stderr, 
		    "%s: error: can't open file: %s\n", progname,instring);
		exit(1);
	    }
            break;
 	case 'q':
	    quiet_flag++;
	    break;
        case '?':
            errflag++;
    }


    if (errflag) {
        fprintf(stderr, "usage: %s [options] < input_file\n\n", progname);
	fprintf(stderr, "    [-t filter_type] \"butterworth\" (default), \"maxflat\", \"chebyshev\"\n");
	fprintf(stderr, "                     \"modcheb\", \"invcheb\", or \"bessel\"\n");
        fprintf(stderr, "    [-n filter_order]             (default = 8)\n");
        fprintf(stderr, "    [-c center_frequency/nyquist] (default = .5)\n");
        fprintf(stderr, "    [-b bandwidth/nyquist]        (default = .1)\n");
	fprintf(stderr, "    [-r ripple_tolerance (dB)]    (default = 3.0103)\n");
        fprintf(stderr, "    [-s Fsb/Fpb] inverse-chebyshev shape factor (default = 2.0)\n");
        fprintf(stderr, "    [-Fs] Normalize frequencies to the sampling rate, rather than Nyquist\n");
	fprintf(stderr, "    [-i dc_value] initialize filter state corresponding to input dc_value\n");
        fprintf(stderr, "    [-x] output magnitude response in autoplot(1) format\n");
        fprintf(stderr, "    [-X] output magnitude and delay plot\n");
        fprintf(stderr, "    [-L] output loss and delay plot\n");
        fprintf(stderr, "    [-P] output loss and delay plot, passband detail\n");
        fprintf(stderr, "    [-I #samples] output impulse response plot\n");
        fprintf(stderr, "    [-S #samples] output step response plot\n");
        fprintf(stderr, "    [-d] dump z-domain biquad values to stdout\n");
        fprintf(stderr, "    [-f file] get z-domain biquad values from file\n");
	fprintf(stderr, "    [-q] do not pass non-numeric lines to the output\n\n");
        fprintf(stderr, "  For example: if the sampling rate is 1 kHz and the desired 3 dB\n");
        fprintf(stderr, "  frequency is 50 Hz, a 20th-order filter would be invoked as:\n\n");
        fprintf(stderr, "    \"%s -b0.1 -n20 < input_file > output_file\"\n\n", progname);
        fprintf(stderr, "  Input format is \"time value\" pairs, one pair per line.\n");
        fprintf(stderr, "  The \"time\" field is ignored, but is echoed to output.\n");
        exit(2);
    }

    if (fnorm == 0.0)
	fnorm = FNYQUIST;  /* default frequency normalization (nyquist) */

    if ((bw<0.0)||(bw>FNYQUIST/fnorm)) {
        fprintf(stderr,
	    "%s: error: bandwidth must be between 0.0 and %.1f\n", progname, FNYQUIST/fnorm);
	exit(1);
    }

    if (filter_type & INV_CHEB && ((sb_ratio<=1.0)||(sb_ratio>=(FNYQUIST/fnorm)/bw))) {
        fprintf(stderr,
	    "%s: error: stopband/passband ratio must be between 1.0 and %.1f/bandwidth\n", progname, FNYQUIST/fnorm);
	exit(1);
    }

    if ((cf<bw/2.0 || cf>(FNYQUIST/fnorm)-bw/2.0) && (filter_type & (BP|BR))) {
        fprintf(stderr,
	    "%s: error: center freq. must be between BW/2 and %.1f-BW/2\n", progname, FNYQUIST/fnorm);
	exit(1);
    }

    if (filter_type & BESSEL && filter_type & (HP|BP|BR)) {
        fprintf(stderr,
	    "%s: warning: Highpass and bandpass/reject transformations do not\npreserve bessel linear phase\n", progname);
    }

    if (fp != NULL) {
	filter_type = UNKNOWN;	    
    }


    /* initialize biquad structures */

    if (fp) {	    /* use external file for z domain info */

	quadlist = getquad(quadlist,fp, &order);

    } else {	    /* calculate our own parameters */

	/* WARNING: Poles and zeros in this section of code are
	 * represented by their RHP reflections, i.e.,
	 * x is really -x, u is really -u.  This is corrected
	 * when x and u are passed to function calls.
	 */

	nquads = 1 + (order-1)/2;
	eps = sqrt(pow(10.0, ripple/10.0) - 1.0);
	mf_radius = pow(eps, -1.0/order);

	Omegax_sq = cos(M_PI*(order-1)/(2*order));
	Omegax_sq *= Omegax_sq;

	if (filter_type & INV_CHEB) {
	    /* adjust bw and eps for inversion
	     * (must account for sb_ratio warping...tricky!) */
	    if (filter_type & LP)
		lambda = prewarp(sb_ratio*bw)/prewarp(bw);
	    else if (filter_type & HP)
		lambda = prewarp(bw)/prewarp(bw/sb_ratio);
	    else if (filter_type & BP)
	      lambda =
		(prewarp(cf+sb_ratio*bw/2.0) - prewarp(cf-sb_ratio*bw/2.0)) /
		(prewarp(cf+bw/2.0) - prewarp(cf-bw/2.0));
	    else if (filter_type & BR)
	      lambda =
		(prewarp(cf+bw/2.0) - prewarp(cf-bw/2.0)) /
		(prewarp(cf+bw/(2.0*sb_ratio)) - prewarp(cf-bw/(2.0*sb_ratio)));
	    else
		lambda = 1.0;

	    eps = 1.0/(eps*ChebPoly(order,lambda));
	}

	if (filter_type & MOD_CHEB && order%2 != 0)
	    filter_type = (filter_type & ~FT_APPROX) | CHEBYSHEV;

	/* Adjust dc gain of even Chebyshevs */
	if (filter_type & CHEBYSHEV && order%2 == 0) {
	    /* set max pband gain to 1.0 */
	    g = 1.0/pow(10.0, ripple/(20.0*nquads));
	} else {
	    g = 1.0;
	}

	for(n=0; n<nquads; n++) { 
	    u = v = 1.0; /* zeros are not used if u != 0.0 */
	    theta=M_PI*(2*(double)n+1)/(2*(double)order);
	    if (filter_type & BUTTERWORTH) {
		x = sin(theta);
		y = cos(theta);
	    } else if (filter_type & MAXFLAT) {
		x = mf_radius*sin(theta);
		y = mf_radius*cos(theta);
	    } else if (filter_type & CHEBYSHEV) {
		x = sin(theta) * sinh(asinh(1.0/eps)/(double) order);
		y = cos(theta) * cosh(asinh(1.0/eps)/(double) order);
	    } else if (filter_type & MOD_CHEB) {
		/* Modified Chebyshev: lowest reflection-zero pair
		 * warped to dc.  See Sedra & Brackett, pp. 144--146.
		 *
		 * s'^2 = (s^2 + Omegax_sq)/(1 - Omegax_sq)
		 */
		double xx, yy;
		x = sin(theta) * sinh(asinh(1.0/eps)/(double) order);
		y = cos(theta) * cosh(asinh(1.0/eps)/(double) order);
		xx = (Omegax_sq + x*x - y*y)/(1-Omegax_sq);
		yy = (2.0*x*y)/(1-Omegax_sq);
		ComplexSqrt(&xx, &yy);
		x = xx;  y = yy;
	    } else if (filter_type & INV_CHEB) {
	        /* See S & B, pp. 147--152
	         * Invert chebyshev poles and add tzeros
	         *
	         * 1/(x+jy) = (x-jy)/(x^2+y^2)
	         */
	        double mag_sq;
		x = sin(theta) * sinh(asinh(1.0/eps)/(double) order);
		y = cos(theta) * cosh(asinh(1.0/eps)/(double) order);
	        mag_sq = x*x+y*y;
	        x =  x*lambda/mag_sq;
	        y = -y*lambda/mag_sq;
	        /* transmission zeros */
	        if ( 2*(n+1) <= order ) {
		    u = 0.0;
		    v = lambda/cos(theta);
	        }
	    } else if (filter_type & BESSEL) {
		if (order > 20) {
		    fprintf(stderr,
		    "%s: error: max bessel filter order = 20\n", progname);
		    exit(1);
		}
		x = bessel[order][2*n];
		y = bessel[order][2*n+1]; 
	    } else {
		fprintf(stderr,
		    "%s: error: filter type not yet implemented\n", progname);
		exit(1);
	    }
	    if ( 2*(n+1) > order ) {	/* on axis, do simple pole */
		quadlist=initpole(filter_type,cf,bw,-x,quadlist);
	    } else {    		/* off axis, do complex pole */
		quadlist=initquad(filter_type,cf,bw,g,-u,v,-x,y,quadlist);
	    }
	}

    }

    initfilter(&dc,quadlist);

    if (dump_flag) {	/* just dump biquad values and quit */
	
	printf("# Dump format is \"c0 c1 c2 / d0 d1 d2\",\n");
	printf("# where the z-domain biquadratic function is:\n");
	printf("#    c0 + c1*z^(-1) + c2*z^(-2)\n");
	printf("#   ----------------------------\n");
	printf("#    d0 + d1*z^(-1) + d2*z^(-2)\n");
	printf("#\n");

	printquad(quadlist);	

    } else if (xfer_plot) {    	/* just do transfer fxn plot */

    	int npoints = 10*order;
    	double F1 = 0.0, F2 = FNYQUIST/fnorm;

    	if (filter_type & (BP|BR))  npoints *= 2;
	if (npoints < 50)  npoints = 50;

	if (xfer_plot & PLOT_PBAND) {
	    npoints = 400;
	    switch (filter_type & FT_BAND) {
	      case LP:  F1 = 0.0;  F2 = bw;
	                break;
	      case HP:  F1 = bw;   F2 = FNYQUIST/fnorm;
	                break;
	      case BP:  F1 = cf-bw/2.0;
	                F2 = cf+bw/2.0;
	                break;
	      case BR:  F1 = 0.0;	/* show just the lower band */
	                F2 = cf-bw/2.0;
	                break;
	      default:  F1 = 0.0;  F2 = FNYQUIST/fnorm;
	    }
	}

	if (xfer_plot & PLOT_MAGN) {
	    initxplot(filter_type,order,cf,bw,"transfer");
	    for (F=F1; F<=F2; F+=bw/(double)npoints) {
		outr=1.0, outi=0.0;
		fevalquad(F,&outr,&outi,quadlist);
		out = sqrt(outi*outi+outr*outr);
		printf("%g %g\n", F, out);
	    }
	} else {
	    initxplot(filter_type,order,cf,bw,"loss");
	    for (F=F1; F<=F2; F+=bw/(double)npoints) {
		outr=1.0, outi=0.0;
		fevalquad(F,&outr,&outi,quadlist);
		out = sqrt(outi*outi+outr*outr);
		out = (out > 1e-20) ? 1.0/out : 1e20;
		printf("%g %g\n", F, out);
	    }
	}

	if (xfer_plot & PLOT_DELAY) {
	    printf("rightyscale 1 delay [samples];; liny\n");
	    for (F=F1; F<=F2; F+=bw/(double)npoints) {
		delay = 0.0;
		devalquad(F,&delay,quadlist);
		printf("%g %g\n", F, delay); 
	    }
	}

    } else if (Isamples) {    	/* output impulse response */

	initIplot(filter_type,order,cf,bw,"impulse");

	out = 1.0;
	filter(&out,quadlist);
	printf("0 %g\n", out);
    	for (n = 1; n < Isamples; n++) {
	    out = 0.0;
	    filter(&out,quadlist);
	    printf("%d %g\n", n, out);
	}

    } else if (Ssamples) {    	/* output step response */

	initIplot(filter_type,order,cf,bw,"step");

    	for (n = 0; n < Ssamples; n++) {
	    out = 1.0;
	    filter(&out,quadlist);
	    printf("%d %g\n", n, out);
	}

    } else { 	    	/* process each x,y pair from standard input */

	while(gets(instring) != NULL) {
    	    scan_OK=sscanf(instring,"%lf %lf",&time, &out);
    	    if(scan_OK == 2) {
	    	filter(&out, quadlist);
	    	printf("%g %g\n", time, out);
            } else if (! quiet_flag) {
		 puts(instring);
	    }
    	}
    }
    exit(0);
}


/* calculates biquad coefficients for complex zero/pole pair */
struct biquad *initquad(filter_type,cf,bw,g,u,v,x,y,p)
int filter_type;	/* type of filter */
double cf,bw;		/* center freq & bandwidth */
double g;	    	/* dc gain of biquad section */
double u,v;	    	/* real and imaginary coordinates of complex zero */
double x,y;	    	/* real and imaginary coordinates of complex pole */
struct biquad *p;   	/* pointer to biquad data structure */
{
    double w0, wL, wH, wn_sq;
    double a, b;
    double x01, y01, x02, y02, v01, v02;
    double f0, alpha;
    struct biquad *q;


    if (filter_type & LP) {

	/* do frequency scale with bilinear warping */
	w0 = prewarp(bw);
	v *= w0; x *= w0; y *= w0;

	wn_sq = x*x + y*y;
    	if ( u != 0.0 )	/* just poles */
    	    p=Ztransform(g*wn_sq, 0.0, 0.0, wn_sq, -2.0*x, 1.0, p);
        else		/* poles and zeros */
    	    p=Ztransform(g*wn_sq, 0.0, g*wn_sq/(v*v), wn_sq, -2.0*x, 1.0, p);

    } else if (filter_type & HP) {

	/* do frequency scale with bilinear warping */
	w0 = prewarp(bw);

	wn_sq = x*x + y*y;
	x = w0*x/wn_sq;
	y = -w0*y/wn_sq;

	wn_sq = u*u + v*v;
	u = w0*u/wn_sq;
	v = -w0*v/wn_sq;

	wn_sq = x*x + y*y;
    	if ( u != 0.0 )	/* just poles */
    	    p=Ztransform(0.0, 0.0, g, wn_sq, -2.0*x, 1.0, p);
        else		/* poles and zeros */
    	    p=Ztransform(g*v*v, 0.0, g, wn_sq, -2.0*x, 1.0, p);

    } else if (filter_type & (BP|BR)) {

	if (filter_type & BR) {
	    /* LP -> HP first */
	    wn_sq = x*x + y*y;
	    x =  x/wn_sq;
	    y = -y/wn_sq;

	    wn_sq = u*u + v*v;
	    u =  u/wn_sq;
	    v = -v/wn_sq;
	}

	/***********************************************************/
	/* compute frequency-warped "analog" specs */

	wL = prewarp(cf-bw/2.0);
	wH = prewarp(cf+bw/2.0);

	/* *geometric* center of warped frequencies */
	w0 = sqrt(wL*wH);

	/* denormalize bandwidth */
	x *= (wH - wL);
	y *= (wH - wL);
	v *= (wH - wL);

	/***********************************************************/
	/* do bandpass transform */

	a = (x*x - y*y - 4.0*w0*w0);
	b = (2.0*x*y);
	ComplexSqrt(&a, &b);

	x01=(x + a)/2.0; y01=(y + b)/2.0;
	x02=(x - a)/2.0; y02=(y - b)/2.0;

	/* only relevant if u == 0 */
	if ( u == 0.0 ) {
	    a = (0.0 - v*v - 4.0*w0*w0);
	    b = 0.0;
	    ComplexSqrt(&a, &b);

	    v01=(v + b)/2.0;
	    v02=(v - b)/2.0;
	}
	/***********************************************************/

	g = sqrt(g);	/* decrease gain per biquad */
	f0 = warp(w0);	/* geometric center frequency */

	wn_sq = x01*x01 + y01*y01;

    	if (filter_type & BP) {
    	    if ( u != 0.0 )	/* just poles */
    	        p=Ztransform(0.0, 1.0, 0.0, wn_sq, -2.0*x01, 1.0, p);
            else		/* poles and zeros */
    	        p=Ztransform(1.0, 0.0, 1.0/(v01*v01), wn_sq, -2.0*x01, 1.0, p);

    	    /* scale gain at center frequency */
    	    q = p->next;
    	    p->next = NULL;
    	    a = 1.0; b = 0.0;
    	    fevalquad(f0, &a, &b, p);
    	    alpha = sqrt(a*a + b*b)/g;
    	    p->b0 /= alpha;
    	    p->b1 /= alpha;
    	    p->b2 /= alpha;
    	    p->next = q;

	} else if (filter_type & BR) {
    	    if ( u != 0.0 )	/* just poles */
    	        p=Ztransform(g*w0*w0, 0.0, g, wn_sq, -2.0*x01, 1.0, p);
            else		/* poles and zeros */
    	        p=Ztransform(g*wn_sq, 0.0, g*wn_sq/(v01*v01), wn_sq, -2.0*x01, 1.0, p);
	}

	wn_sq = x02*x02 + y02*y02;

    	if (filter_type & BP) {
    	    if ( u != 0.0 )	/* just poles */
    	        p=Ztransform(0.0, 1.0, 0.0, wn_sq, -2.0*x02, 1.0, p);
            else		/* poles and zeros */
    	        p=Ztransform(1.0, 0.0, 1.0/(v02*v02), wn_sq, -2.0*x02, 1.0, p);

    	    /* scale gain at center frequency */
    	    q = p->next;
    	    p->next = NULL;
    	    a = 1.0; b = 0.0;
    	    fevalquad(f0, &a, &b, p);
    	    alpha = sqrt(a*a + b*b)/g;
    	    p->b0 /= alpha;
    	    p->b1 /= alpha;
    	    p->b2 /= alpha;
    	    p->next = q;

	} else if (filter_type & BR) {
    	    if ( u != 0.0 )	/* just poles */
    	        p=Ztransform(g*w0*w0, 0.0, g, wn_sq, -2.0*x02, 1.0, p);
            else		/* poles and zeros */
    	        p=Ztransform(g*wn_sq, 0.0, g*wn_sq/(v02*v02), wn_sq, -2.0*x02, 1.0, p);
	}

    }

    return(p);
}

/* calculates biquad coefficients for simple pole on x axis */
struct biquad *initpole(filter_type,cf,bw,x,p)
int filter_type;	/* type of filter */
double cf,bw;		/* center freq & bandwidth */
double x;	    	/* x coordinate of real pole */
struct biquad *p;   	/* pointer to biquad data structure */
{
    double w0, wL, wH, wn_sq;
    double a, b, y, x01, y01, x02, y02;
    double f0, alpha;
    struct biquad *q;


    if (filter_type & LP) {

	/* do frequency scale with bilinear warping */
	w0 = prewarp(bw);
	x *= w0; 

	p=Ztransform(1.0, 0.0, 0.0, 1.0, -1.0/x, 0.0, p);

    } else if (filter_type & HP) {

	/* do frequency scale with bilinear warping */
	w0 = prewarp(bw);
	x = w0/x; 

	p=Ztransform(0.0, -1.0/x, 0.0, 1.0, -1.0/x, 0.0, p);

    } else if (filter_type & (BP|BR)) { 

	if (filter_type & BR) {
	    /* LP -> HP first */
	    x = 1.0/x;
	}

	/***********************************************************/
	/* compute frequency-warped "analog" specs */

	wL = prewarp(cf-bw/2.0);
	wH = prewarp(cf+bw/2.0);

	/* *geometric* center of warped frequencies */
	w0 = sqrt(wL*wH);

	/* denormalize bandwidth */
	x *= (wH - wL);
	y = 0.0;

	/***********************************************************/
	/* do bandpass transform */

	a = x*x - y*y - 4.0*w0*w0;	/* x^2 - 4 w0^2 */
	b = 2.0*x*y;			/* == 0 */

	if (a < 0.0) {
	    /* conjugate roots */
	    ComplexSqrt(&a, &b);
	    x02 =   x01 = x/2.0;	/* (x +- a)/2, but a==0 */
	    y02 = -(y01 = b/2.0);	/* (y +- b)/2, but y==0 */
	} else {
	    /* real roots */
	    a = sqrt(a);
	    x01 = (x + a)/2.0;
	    x02 = (x - a)/2.0;
	    y02 = y01 = 0.0;
	}

	/***********************************************************/

	f0 = warp(w0);	/* geometric center frequency */

	wn_sq = x01*x02 - y01*y02;

    	if (filter_type & BP) {
    	    p=Ztransform(0.0, 1.0, 0.0, wn_sq, -(x01+x02), 1.0, p);

    	    /* scale gain at center frequency */
    	    q = p->next;
    	    p->next = NULL;
    	    a = 1.0; b = 0.0;
    	    fevalquad(f0, &a, &b, p);
    	    alpha = sqrt(a*a + b*b);
    	    p->b0 /= alpha;
    	    p->b1 /= alpha;
    	    p->b2 /= alpha;
    	    p->next = q;

	} else if (filter_type & BR) {
    	    p=Ztransform(w0*w0, 0.0, 1.0, wn_sq, -(x01+x02), 1.0, p);
	}
    }
    return(p);
}

double filter(vin,p)
double *vin;
struct biquad *p;
{
    double v0, vout;

    /*
    ** compute recursive digital filter given coefficients
    ** a1,a2,b0,b1,b2, with intermediate storage values v1, v2,
    ** all stored in "biquad" structure.
    **
    **                         (v0)   
    **  Vin -->(+)-----------.--->[b0]---->(+)-----> Vout
    **          ^            |              ^
    **          |          [1/z]            |
    **          |            | (v1)         |
    **          <----[a1]<---.--->[b1]------>
    **          |            |              |
    **          |          [1/z]            |
    **          |            | (v2)         |
    **          <----[a2]<---.--->[b2]------>
    **
    **
    **  Vout       b0 + b1*z^(-1) + b2*z^(-2)
    ** ------  =  ----------------------------
    **  Vin         1 - a1*z^(-1) - a2*z^(-2)
    */

    if (p != NULL) {
    	filter(vin,p->next);

	v0 = *vin + p->a1*p->v1 + p->a2*p->v2;
	vout = p->b0*v0 + p->b1*p->v1 + p->b2*p->v2;

	/* shift samples in time */
	p->v2 = p->v1;
	p->v1 = v0;

	*vin = vout;
    } 
}

double initfilter(vdc,p)
double *vdc;
struct biquad *p;
{
    double v0, vout;

    /*
    ** Initialize digital filter state as if it has "seen"
    ** a dc input of vdc since -infinity.  Note that in steady
    ** state (dc), v2 = v1 = v0.  See filter() for diagram of
    ** filter structure.
    **
    ** sdw, 7/30/93
    **
    */

    if (p != NULL) {
    	initfilter(vdc,p->next);

	v0 = *vdc / (1.0 - p->a1 - p->a2);
	vout = v0 * (p->b0 + p->b1 + p->b2);

	/* copy states */
	p->v2 = v0;
	p->v1 = v0;

	*vdc = vout;
    } 
}

struct biquad *getquad(p,fp,order)  /* get biquad values from file */
struct biquad *p;
FILE *fp;
int *order;
{
    struct biquad *q;  	/* pointer to biquad data structure */
    double c0,c1,c2,d0,d1,d2;
    int n = 0;

    while (fgets(instring,MAXBUF,fp) != NULL) {
    	scan_OK=sscanf(instring,"%lf %lf %lf / %lf %lf %lf",
	    	&c0, &c1, &c2, &d0, &d1, &d2);

    	if(scan_OK==6) {

	    /* update filter order, checking if bq is 1st order */
	    if (c2 == 0.0 && d2 == 0.0)
		n += 1;
	    else
		n += 2;

  	    q = qalloc();        /* make a new node */
	    q->next = p;         /* stick old list after this one */
	    p = q;		 /* set p to new head */

	    q->v1 = 0.0;
	    q->v2 = 0.0;

	    if (d0 == 0.0) {
        	fprintf(stderr, 
		    "%s: error: d0 coefficient must be non-zero\n", progname);
		exit(1);
	    }

	    q->b0 = c0/d0;
	    q->b1 = c1/d0;
	    q->b2 = c2/d0;

	 /* q->a0 = d0/d0;  a0 is not in structure, just assumed */

	    q->a1 = -d1/d0;
	    q->a2 = -d2/d0;

        } 
    }

    *order = n;
    return(p);
}

printquad(p)	    /* print out values of biquad structure */
struct biquad *p;
{
    if (p != NULL) {
    	printquad(p->next);
	printf("%-.10g %-.10g %-.10g / %-.10g %-.10g %-.10g\n",
	    p->b0, p->b1, p->b2, 1.0, -p->a1, -p->a2);
    }
}

initxplot(filter_type,order,cf,bw,ptype)
int filter_type;
int order;
double cf,bw;
char *ptype;	/* "transfer" or "loss" */
{
    char *name;

    switch (filter_type & FT_APPROX) {
	case BUTTERWORTH:	name = "Butterworth";
				break;
	case MAXFLAT:		name = "Maximally Flat";
				break;
	case CHEBYSHEV:		name = "Chebyshev";
				break;
	case MOD_CHEB:		name = "Modified Chebyshev";
				break;
	case INV_CHEB:		name = "Inverse Chebyshev";
				break;
	case BESSEL:		name = "Bessel";
				break;
	default:		name = "";
    }

    switch (filter_type & FT_BAND) {
	case LP:
    	    printf("title %s lowpass %s function\n",name,ptype);
    	    printf("title BW = %g, Order = %d\n",bw,order);
	    break;
	case HP:
    	    printf("title %s highpass %s function\n",name,ptype);
    	    printf("title BW = %g, Order = %d\n",bw,order);
	    break;
	case BP:
    	    printf("title %s bandpass %s function\n",name,ptype);
    	    printf("title Fcenter = %g, BW = %g, Order = %d\n",cf,bw,order);
	    break;
	case BR:
    	    printf("title %s bandreject %s function\n",name,ptype);
    	    printf("title Fcenter = %g, BW = %g, Order = %d\n",cf,bw,order);
	    break;
	default:
    	    printf("title filter %s function\n",ptype);
    	    printf("title Order = %d\n",order);
    }

    printf("style presentation\n");
    if (fnorm == FSAMPLE)
	printf("xscale 1 frequency [normalized to Fs]\n");
    else if (fnorm == FNYQUIST)
	printf("xscale 1 frequency [normalized to Fs/2]\n");
    else
	printf("xscale 1 frequency [normalized to %.2f*Fs]\n", fnorm);
    if ( ptype && (*ptype == 'l' || *ptype == 'L') )
	printf("yscale 1 loss [dB];; dby\n");
    else
	printf("yscale 1 magnitude [dB];; dby\n");
}

initIplot(filter_type,order,cf,bw,ptype)
int filter_type;
int order;
double cf,bw;
char *ptype;	/* "impulse" or "step" */
{
    char *name;

    switch (filter_type & FT_APPROX) {
	case BUTTERWORTH:	name = "Butterworth";
				break;
	case MAXFLAT:		name = "Maximally Flat";
				break;
	case CHEBYSHEV:		name = "Chebyshev";
				break;
	case MOD_CHEB:		name = "Modified Chebyshev";
				break;
	case INV_CHEB:		name = "Inverse Chebyshev";
				break;
	case BESSEL:		name = "Bessel";
				break;
	default:		name = "";
    }

    switch (filter_type & FT_BAND) {
	case LP:
    	    printf("title %s lowpass %s response\n",name,ptype);
    	    printf("title BW = %g, Order = %d\n",bw,order);
	    break;
	case HP:
    	    printf("title %s highpass %s response\n",name,ptype);
    	    printf("title BW = %g, Order = %d\n",bw,order);
	    break;
	case BP:
    	    printf("title %s bandpass %s response\n",name,ptype);
    	    printf("title Fcenter = %g, BW = %g, Order = %d\n",cf,bw,order);
	    break;
	case BR:
    	    printf("title %s bandreject %s response\n",name,ptype);
    	    printf("title Fcenter = %g, BW = %g, Order = %d\n",cf,bw,order);
	    break;
	default:
    	    printf("title filter %s response\n",ptype);
    	    printf("title Order = %d\n",order);
    }

    printf("style working\n");
    printf("xscale 1 samples\n");
    printf("yscale 1 amplitude\n");
    printf("symbol circle;; symbolsize 0.5\n");
}

double fevalquad(F,realp,imagp,p) /* realp, imagp are _inputs_ and outputs */
double F;
double *realp, *imagp;
struct biquad *p;
{
    double real, imag, wT;
    double sin1, sin2, cos1, cos2, mag_sq;
    double nreal = 0.0; double nimag = 0.0;
    double dreal = 1.0; double dimag = 0.0;

    if (p != NULL) {
    	fevalquad(F,realp,imagp,p->next);
	
	wT = 2.0*M_PI*(F*fnorm);

#if 0
	nreal += p->b0*cos(-0.0*wT);	nimag += p->b0*sin(-0.0*wT);			
	nreal += p->b1*cos(-1.0*wT);	nimag += p->b1*sin(-1.0*wT);    
	nreal += p->b2*cos(-2.0*wT);	nimag += p->b2*sin(-2.0*wT);    

	dreal -= p->a1*cos(-1.0*wT);	dimag -= p->a1*sin(-1.0*wT);
	dreal -= p->a2*cos(-2.0*wT);	dimag -= p->a2*sin(-2.0*wT);   
#else
	/* a little speedup */
	sin1 = -sin(wT);  sin2 = -sin(2.0*wT);
	cos1 =  cos(wT);  cos2 =  cos(2.0*wT);

	nreal += p->b0;
	nreal += p->b1*cos1;	nimag += p->b1*sin1;    
	nreal += p->b2*cos2;	nimag += p->b2*sin2;    

	dreal -= p->a1*cos1;	dimag -= p->a1*sin1;
	dreal -= p->a2*cos2;	dimag -= p->a2*sin2;   
#endif

	/* now divide complex numerator by denominator */
	mag_sq = dimag*dimag + dreal*dreal + 1e-20;
	real = (nreal*dreal + nimag*dimag)/mag_sq;
	imag = (nimag*dreal - nreal*dimag)/mag_sq;

	/* and multiply input pointers by the result, (reusing dreal,dimag) */
	dreal = (*realp) * real - (*imagp) * imag;
	dimag = (*realp) * imag + (*imagp) * real;

	/* now assign back to input pointers */
	*realp = dreal;
	*imagp = dimag;
    }
}

double devalquad(F,Tau,p)  /* Tau is _input_ and output, normalized by T */
double F;
double *Tau;
struct biquad *p;
{
    double wT;
    double c0, c1, c2;
    double sin1, sin2, cos1, cos2;
    double N, D, dN_dw, dD_dw, SS;

    /* Example delay computation:
     *
     * Delay of zero polynomial c0 + c1*z^(-1) + c2*z^(-2)
     * where z = exp(jwT) and wT = pi*f/fN = pi*F
     *
     * Phase is:  phi(w) = atan(N/D)
     *     where  N = -c1*sin(wT) - c2*sin(2wT)
     *       and  D = c0 + c1*cos(wT) + c2*cos(2wT)
     *
     * Delay is:  tau(w) = -dphi/dw
     *                   = (N*dD/dw - D*dN/dw) / (N^2+D^2)
     *
     *      and:  dN/dw = -c1*T*cos(wT) - 2*c2*T*cos(2wT)
     *            dD/dw = -c1*T*sin(wT) - 2*c2*T*sin(2wT)
     */

    if (p != NULL) {
    	devalquad(F,Tau,p->next);
	
	wT = 2.0*M_PI*(F*fnorm);

	/* avoid recomputations */
	sin1 = sin(wT);  sin2 = sin(2.0*wT);
	cos1 = cos(wT);  cos2 = cos(2.0*wT);

	/* delay of zeros */
	c0 = p->b0; c1 = p->b1; c2 = p->b2;
	N = 0.0 - c1*sin1 - c2*sin2;
	D =  c0 + c1*cos1 + c2*cos2;
	dN_dw = 0.0 - c1*T*cos1 - 2.0*c2*T*cos2;
	dD_dw = 0.0 - c1*T*sin1 - 2.0*c2*T*sin2;

	if ( (SS = N*N + D*D) != 0.0 )
	    *Tau += (N*dD_dw - D*dN_dw) / (T*SS);

	/* delay of poles */
	c0 = 1.0; c1 = -p->a1; c2 = -p->a2;
	N = 0.0 - c1*sin1 - c2*sin2;
	D =  c0 + c1*cos1 + c2*cos2;
	dN_dw = 0.0 - c1*T*cos1 - 2.0*c2*T*cos2;
	dD_dw = 0.0 - c1*T*sin1 - 2.0*c2*T*sin2;

	if ( (SS = N*N + D*D) != 0.0 )
	    *Tau += (D*dN_dw - N*dD_dw) / (T*SS);
    }
}

double prewarp(F)	/* bilinear prewarp: digital F to analog w */
double F;
{
    return( (2.0/T)*tan(M_PI*F*fnorm) );
}

double warp(w)		/* bilinear warp: analog w to digital F */
double w;
{
    return ( (1.0/(M_PI*fnorm))*atan(w*T/2.0) );
}

struct biquad *Ztransform(a0,a1,a2,b0,b1,b2,p) /* general 2nd-order biquad */
double a0,a1,a2,b0,b1,b2;
struct biquad *p;   	/* pointer to biquad data structure */
{
    double c0,c1,c2,d0,d1,d2;
    struct biquad *q;  	/* pointer to biquad data structure */

    /* these equations are computed by taking an equation
    ** with both quadratic numerator and denominator:
    **
    **      Vout     [a0 + a1*s + a2*s^2] 
    **     ------ = ----------------------
    **      Vin      [b0 + b1*s + b2*s^2]
    **
    ** and using the bilinear transform to convert to z-domain:
    **
    **                2 * (z-1)
    **       s    =  -----------
    **                T * (z+1)
    **
    ** Where z^(-n) is a delay of n time steps, and T = sample timestep.
    **
    ** sdw, 8/12/93: This is modified when a2==b2==0.
    */

    if (a2 != 0.0 || b2 != 0.0) {
	c0 = 4.0*a2 + 2.0*a1*T + a0*T*T;
	c1 = 2.0*a0*T*T - 8.0*a2;       
	c2 = 4.0*a2 - 2.0*a1*T + a0*T*T;

	d0 = 4.0*b2 + 2.0*b1*T + b0*T*T;
	d1 = 2.0*b0*T*T - 8.0*b2;       
	d2 = 4.0*b2 - 2.0*b1*T + b0*T*T;
    } else if (a1 != 0.0 || b1 != 0.0) {
	c0 = a0*T + 2.0*a1;
	c1 = a0*T - 2.0*a1;
	c2 = 0.0;

	d0 = b0*T + 2.0*b1;
	d1 = b0*T - 2.0*b1;
	d2 = 0.0;
    } else {
	c0 = a0; c1 = 0.0, c2 = 0.0;
	d0 = b0; d1 = 0.0, d2 = 0.0;
    }

    q = qalloc();	 /* make a new node */
    q->next = p;	 /* stick old list after this one */

    /* Note: biquad a's and b's differ from those above!
     * see filter()
     */
    q->b0 = c0/d0;
    q->b1 = c1/d0;
    q->b2 = c2/d0;
    q->a1 = -d1/d0;
    q->a2 = -d2/d0;
    q->v1 = 0.0;
    q->v2 = 0.0;
    
    return(q);
}

struct biquad *qalloc()
{
    return((struct biquad *) malloc((unsigned) sizeof(struct biquad)));
}

ComplexSqrt(x,y)    /* complex square root */
double *x, *y;
{
    double rho,theta;

    rho = sqrt(sqrt((*x) * (*x) + (*y) * (*y)));
    theta = atan2((*y),(*x))/2.0;
    *x = rho*cos(theta);
    *y = rho*sin(theta);
}

double asinh(x)    /* arcsinh for real values */
double x;
{
    return(log(x + sqrt(x*x + 1.0)));
}

double ChebPoly(N, x)	/* value of Nth-order Chebyshev polynomial at x */
int N;
double x;
{
    if (N < 0) {
	fprintf(stderr, 
	    "%s: error in ChebPoly(): N must be greater than 0\n", progname);
	return(0.0);
    } else if (N == 0) {
        return(1.0);
    } else if (N == 1) {
        return(x);
    } else {
        return( 2*x*ChebPoly(N-1, x) - ChebPoly(N-2, x) );
    }
}
