#ifndef __MATHUTIL_H__
#define __MATHUTIL_H__

#ifdef __cplusplus
extern "C" {
#endif

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0F ? 0.0F : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static float cubarg;
#define CUB(a) ((cubarg=(a)) == 0.0F ? 0.0F : cubarg*cubarg*cubarg)

static double dcubarg;
#define DCUB(a) ((dcubarg=(a)) == 0.0 ? 0.0 : dcubarg*dcubarg*dcubarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int ipararg;
#define PAR(a) (ipararg=(a), (ipararg)/ 2 << 1 != (ipararg) ? -1 : 1 )

/***********************************************************/
/*                                                         */
/*   Angular momentum coeffcients calcualtion              */
/*                                                         */
/*    Previously used at PDP-11 computer, modified         */
/*   to work at IBM - AT                                   */
/***********************************************************/
/* translated by f2c (version of 21 October 1993  13:46:10)*/
/*   end edited by Sergey Panov  Jan 4 1994                */
/***********************************************************/

double c_(double *a, double *b, double *g, double *d, 
	double *e, double *f);
double cg_(double *a, double *d);
double del_(double *a, double *b, double *g);
double flog_(double *a);

/*      DOUBLE PRECISION RACAH COEFFICIENTS */
double w_(	double *a, double *b, double *c,
			double *d, double *e, double *f );

/*      DOUBLE PRECISION 3-J SYSMBOL */
double TJS(	double *a, double *b, double *g,
			double *d, double *e, double *f );

/*      DOUBLE PRECISION 6-J SYSMBOL */
/*                        A B E */
/*                        D C F */
double SJS( double a, double b, double e,
			double d, double c, double f );

/*      DOUBLE PRECISION 9-J SYSMBOL */
/*                  A B C */
/*                  D E F */
/*                  H P Q */

double WNJ(	double *a, double *b, double *c,
			double *d, double *e, double *f,
			double *g, double *h, double *p );

/* ******************************************************************** */
/* FUNCTION THAT CALCULATES CLEBSCH-GORDAN COEFFICIENTS. THE FORMULAE   */
/* CAN BE FOUND ON PAGE 57 OF RICHARD N. ZARE'S `ANGULAR MOMENTUM'      */
/* ******************************************************************** */
double CG1(	double j1, double j2, double j3,
			double p1, double p2, double p3 );


#ifdef __cplusplus
}
#endif

#endif //__MATHUTIL_H__
