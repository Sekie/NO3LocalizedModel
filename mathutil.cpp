#include "msdef.h"

#include <math.h>
#include "mathutil.h"

/***********************************************************/
/*                                                         */
/*   Angular momentum coeffcients calcualtion              */
/*                                                         */
/*    Previously used at PDP-11 computer, modified         */
/*   to work at IBM - AT                                   */
/*                                                         */
/*                                                         */
/*   Note:                                                 */
/*         Since most subroutines in this pakage are       */
/*      DOUBLE PRECISION, when calling these routines      */
/*      directly by substituting numbers for the variables */
/*      the number MUST BE EXPRESS IN DOUBLE PRECISION.    */
/*      Otherwise, the results will NOT be correct !       */
/*  EXAMPLE: If one want calculate  3J-symbol of           */
/*           J1 = 2, m1 = 1, J2 = 2  m2 = 2, J3 = 3 m3=-3  */
/*      One should call the function TJS as                */
/*     RESULT = TJS(2.0d0,2.0d0,3.0d0,1.0d0,2.0d0,3.0d0)   */
/*                                                         */
/*                   Xianming Liu                          */
/*                                 Jun-12-87               */
/***********************************************************/
/* translated by f2c (version of 21 October 1993  13:46:10)*/
/*   end edited by Sergey Panov  Jan 4 1994                */

double c_(double *a, double *b, double *g, double *d, 
	double *e, double *f)
{
    /* Format strings */

    /* System generated locals */
    int i__1;
    double r__1, r__2;
    double ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    extern double flog_(double *);
    static int i, j, k, l, m, n;
    static double q, r, s;
    extern double cg_(double *, double *), del_(double *, 
	    double *, double *);

    ret_val = 0.;
    i = (int) (*a + *b + *g);
//    VERIFY((double) i == *a + *b + *g);
//    VERIFY(!(*a < 0. || *b < 0. || ret_val < 0.));
    
    i = (int) (*a + *d);
    j = (int) (*b + *e);
    k = (int) (*g + *f);
//    VERIFY(!((double) i != *a + *d || (double) j != *b + *e || (double)
//          k != *g + *f) );
    i = (int) (*a * 2.);
    j = (int) (*b * 2.);
    k = (int) (*g * 2.);
    l = (int) (*d * 2.);
    m = (int) (*e * 2.);
    n = (int) (*f * 2.);
//    VERIFY( !( (double) i != *a * 2. || (double) j != *b * 2.
//    		|| (double) k != *g * 2. || (double) l != *d * 2. 
//	     	|| (double) m != *e * 2. || (double) n != *f * 2.) );
	if (*d + *e - *f != 0.) return ret_val;
    if (fabs(*d) > *a || fabs(*e) > *b || fabs(*f) > *g) return ret_val;
    if (*a > *b + *g || *b > *a + *g || *g > *a + *b) return ret_val;
    if (*a < (d__1 = *b - *g, fabs(d__1)) || *b < (d__2 = *a - *g, fabs(d__2)) 
	    || *g < (d__3 = *a - *b, fabs(d__3))) return ret_val;
/* Computing MAX */
    r__1 = 0.;
    r__2 = *b - *g - *d;
    r__1 = (r__1 >= r__2) ? r__1 : r__2;
    r__2 = *a - *g + *e;
    l = (int)( (r__1 >= r__2) ? r__1 : r__2 );
    
/* Computing MIN */
    r__1 = *a - *d;
    r__2 = *b + *e;
    r__1 = (r__1 <= r__2) ? r__1 : r__2;
    r__2 = *a + *b - *g;
    m = (int)((r__1 <= r__2) ? r__1 : r__2 );
    
    d__1 = fabs(*d);
    d__2 = fabs(*e);
    d__3 = fabs(*f);
    r = log(*g * 2. + 1.) * .5 + cg_(a, &d__1) + cg_(b, &d__2) + cg_(g, &d__3)
	     + del_(a, b, g);
    d__1 = *a - *d - l;
    d__2 = *g - *b + *d + l;
    d__3 = *b + *e - l;
    d__4 = *g - *a - *e + l;
    d__5 = (double) l;
    d__6 = *a + *b - *g - l;
    r = r - flog_(&d__1) - flog_(&d__2) - flog_(&d__3) - flog_(&d__4) - flog_(
	    &d__5) - flog_(&d__6);
    n = 1;
    if (l / 2 << 1 != l) {
	  n = -1;
    }
    if (l == m) {
      ret_val = n * exp(r);
      return ret_val;
    }
    ret_val = 1.;
    q = 1.;
    i__1 = m;
    for (i = l + 1; i <= i__1; ++i) {
	s = (double) i - 1.;
	q = -q * (*a - *d - s) * (*b + *e - s) * (*a + *b - *g - s) / (*g - *
		b + *d + 1. + s) / (s + 1.) / (*g - *a - *e + (
		float)1. + s);
	ret_val += q;
    }
    if (ret_val == 0.) {
	return ret_val;
    }
    if (ret_val < 0.) {
	n = -n;
    }
    ret_val = n * exp(r + log((fabs(ret_val))));
    return ret_val;
} /* c_ */


double cg_(double *a, double *d)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i;
    static double s;

    ret_val = 0.;
    if (*a == 0.) return ret_val;
    if (*d != *a) {
    	i__1 = (int) (*a - *d);
    	for (i = 1; i <= i__1; ++i) {
			s = (double) i;
			ret_val += log(s);
    	}
    	if (*d == 0.) return ret_val;
    }
    for (i = (int) (*a - *d + 1); i <= i__1; ++i) {
		s = (double) i;
		ret_val += log(s) * .5;
    }
    return ret_val;
} /* cg_ */


double del_(double *a, double *b, double *g)
{
    /* System generated locals */
    int i__1;
    double r__1, r__2;
    double ret_val;

    /* Local variables */
    static int i, l, m, n;
    static double s;

    ret_val = 0.;
/* Computing MIN */
    r__1 = *a + *b - *g;
    r__2 = *a + *g - *b;
    r__1 = ( r__1 <= r__2 ) ? r__1 : r__2 ;
    r__2 = *b + *g - *a;
    l = (int) ((r__1 <= r__2) ? r__1 : r__2);
    
/* Computing MAX */
    r__1 = *a + *b - *g;
    r__2 = *a + *g - *b;
    r__1 = ( r__1 >= r__2 ) ? r__1 : r__2 ;
    r__2 = *b + *g - *a;
    n = (int) ((r__1 >= r__2) ? r__1 : r__2);
    
    m = (int) (*a + *b + *g - l - n);
    
    if (l != 0) {
    	for (i = 1; i <= l; ++i) {
			s = (double) i;
			ret_val += log(s);
    	}
    }
    
    if (m != l) {
    	for (i = l + 1; i <= m; ++i) {
			s = (double) i;
			ret_val += log(s) * .5;
    	}
    }
    
    i__1 = (int) (*a + *b + *g + 1.);
    for (i = n + 1; i <= i__1; ++i) {
		s = (double) i;
		ret_val -= log(s) * .5;
    }
    
    return ret_val;
} /* del_ */


double flog_(double *a)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i;
    static double s;

    ret_val = 0.;
//    VERIFY(*a >= 0.);
    
    if (*a == 0.) return ret_val;
    i__1 = (int) (*a);
    for (i = 1; i <= i__1; ++i) {
		s = (double) i;
		ret_val += log(s);
    }
    return ret_val;
} /* flog_ */


/*      DOUBLE PRECISION RACAH COEFFICIENTS */

double w_(double *a, double *b, double *c, double *d, 
	double *e, double *f)
{
    /* System generated locals */
    int i__1;
    double r__1, r__2;
    double ret_val, d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    extern double flog_(double *);
    static int i, j, k, l, m, n;
    static double q, r, s;
    extern double del_(double *, double *, double *);


    ret_val = 0.;
//    VERIFY(! (*a < 0. || *b < 0. || *c < 0. || *d < 0. || *e < 0. || *f < 0.));
    i = (int) (*a + *b + *e);
    j = (int) (*a + *c + *f);
    k = (int) (*b + *d + *f);
    l = (int) (*c + *d + *e);
//    VERIFY(! ((double) i != *a + *b + *e || (double) j != *a + *c + *f ||
//     			(double) k != *b + *d + *f || (double) l != *c + *d + *e));
    i = (int) (*a * 2.);
    j = (int) (*b * 2.);
    k = (int) (*c * 2.);
    l = (int) (*d * 2.);
    m = (int) (*e * 2.);
    n = (int) (*f * 2.);
//    VERIFY(! ((double) i != *a * 2. || (double) j != *b * 2. 
//	    || (double) k != *c * 2. || (double) l != *d * (
//	    float)2. || (double) m != *e * 2. || (double) n != 
//	    *f * 2.));
    if (*a > *b + *e || *a > *c + *f || *b > *d + *f || *c > *d + *e) {
		return ret_val;
    } else if (*b > *a + *e || *c > *a + *f || *d > *b + *f || *d > *c + *e) {
		return ret_val;
    } else if (*e > *a + *b || *f > *a + *c || *f > *b + *d || *e > *c + *d) {
		return ret_val;
    } else if (*a < (d__1 = *b - *e, fabs(d__1)) || *a < (d__2 = *c - *f, fabs(
	    d__2)) || *b < (d__3 = *d - *f, fabs(d__3)) || *c < (d__4 = *d - *
	    e, fabs(d__4))) {
		return ret_val;
    } else if (*b < (d__1 = *a - *e, fabs(d__1)) 
									|| *c < (d__2 = *a - *f, fabs(d__2))
	    							|| *d < (d__3 = *b - *f, fabs(d__3))
	    							|| *d < (d__4 = *c - *e, fabs(d__4))) {
		return ret_val;
    } else if (*e < (d__1 = *a - *b, fabs(d__1))
									|| *f < (d__2 = *a - *c, fabs(d__2))
									|| *f < (d__3 = *b - *d, fabs(d__3))
									|| *e < (d__4 = *c - *d, fabs(d__4))) {
		return ret_val;
    }
    
    r = del_(a, b, e) + del_(a, c, f) + del_(b, d, f) + del_(c, d, e);
/* Computing MAX */
    r__1 = 0.,	r__2 = *a + *d - *e - *f;
    r__1 = (r__1 >= r__2) ? r__1 : r__2,	r__2 =*b + *c - *e - *f;
    l = (int)((r__1 >= r__2) ? r__1 : r__2);
    
/* Computing MIN */
    r__1 = *a + *b - *e,	r__2 = *c + *d - *e;
    r__1 = (r__1 <= r__2) ? r__1 : r__2,	r__2 = *a + *c - *f;
    r__1 = (r__1 <= r__2) ? r__1 : r__2,	r__2 = *b + *d - *f;
    m = (int)((r__1 <= r__2) ? r__1 : r__2);
    
    d__1 = *a + *b + *c + *d + 1 - l;
    d__2 = (double) l;
    d__3 = *e + *f - *a - *d + l;
    d__4 = *e + *f - *b - *c + l;
    d__5 = *a + *b - *e - l;
    d__6 = *c + *d - *e - l;
    d__7 = *a + *c - *f - l;
    d__8 = *b + *d - *f - l;
    r = r + flog_(&d__1) - flog_(&d__2) - flog_(&d__3) - flog_(&d__4) - flog_(
	    &d__5) - flog_(&d__6) - flog_(&d__7) - flog_(&d__8);
    n = 1;
   	if (l / 2 << 1 != l) n = -1;
    if (m == l) {
    	ret_val = n * exp(r);
    	return ret_val;
    }
    
    ret_val = 1.;
    q = 1.;
    i__1 = m;
    
    for (i = l + 1; i <= i__1; ++i) {
		s = (double) i - 1.;
		q = -q / (*a + *b + *c + *d + 1 - s) / (s + 1) / 
		(*e + *f - *a - *d + s + 1) / (*e + *f - *b - *c + s + 1) *
		(*a + *b - *e - s) * (*c + *d - *e - s) * (*a + *c - *f - s) *
		(*b + *d - *f - s);
		ret_val += q;
    }
    if (ret_val == 0.) 	return ret_val;
    if (ret_val < 0.) n = -n;
    
    ret_val = n * exp(r + log((fabs(ret_val))));
    
    return ret_val;
} /* w_ */


/*      DOUBLE PRECISION 3-J SYSMBOL */

double TJS(double *a, double *b, double *g, double *d, 
	double *e, double *f)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    extern double c_(double *, double *, double *, double*,
     double *, double *);
    static int i, n;

    n = 1;
    i = (int) (*a - *b - *f);
    if (i / 2 << 1 != i)n = -1;
    d__1 = -(*f);
    ret_val = n / sqrt(*g * 2. + 1.) * c_(a, b, g, d, e, &d__1);
    return ret_val;
} /* tjs_ */


/*      DOUBLE PRECISION 6-J SYSMBOL */
/*                        A B E */
/*                        D C F */

double SJS(	double a, double b, double e,
			double d, double c, double f )
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static int i, n;
    extern double w_(double *, double *, double *, double 
	    *, double *, double *);

    n = 1;
    i = (int) (a + b + c + d);
    if (i / 2 << 1 != i) n = -1;
    ret_val = n * w_(&a, &b, &c, &d, &e, &f);
    return ret_val;
} /* SJS */


double SWNJ(double *a, double *b, double *c, double *d, 
	double *e, double *f, double *g, double *h, 
	double *p)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    extern double WNJ(double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *);

    ret_val = WNJ(a, b, c, d, e, f, g, h, p);
    return ret_val;
} /* SWNJ */


/*      DOUBLE PRECISION 9-J SYSMBOL */

/*                  A B C */
/*                  D E F */
/*                  H P Q */

double WNJ(double *a, double *b, double *c, double *d, 
	double *e, double *f, double *g, double *h, 
	double *p)
{
    /* System generated locals */
    int i__1;
    double ret_val, d__1, d__2, d__3, d__4, d__5;

    /* Local variables */
    static int i, j, k, l, m, n;
    static double q;
    extern double w_(double *, double *, double *, double 
	    *, double *, double *);
    static double ak, al, am;
    static int ii, jj, kk;

    ret_val = 0.;
//    VERIFY(! (*a < 0. || *b < 0. || *c < 0. || *d < 0. || *e < 0. || *f < 0. 
//    				  || *g < 0. || *h < 0. || *p < 0.));
	    
    i = (int) (*a * 2.);
    j = (int) (*b * 2.);
    k = (int) (*c * 2.);
    l = (int) (*d * 2.);
    m = (int) (*e * 2.);
    n = (int) (*f * 2.);
    ii = (int) (*g * 2.);
    jj = (int) (*h * 2.);
    kk = (int) (*p * 2.);
//   VERIFY( !((double) i != *a * 2. ||  (double)j != *b * 2. 
//	    || (double) k != *c * 2. || (double) l != *d * 2.
//	    || (double) m != *e * 2. || (double) n != *f * 2.
//	    || (double) ii != *g * 2. || (double) jj != *h * 2.
//	    || (double) kk != *p * 2.) );
	    
    i = (int) ((*a + *p) * 2.);
    j = (int) ((*b + *f) * 2.);
    k = (int) ((*h + *d) * 2.);
    n = 0;
    if (i / 2 << 1 != i) ++n;
    if (j / 2 << 1 != j) ++n;
    if (k / 2 << 1 != k) ++n;
//    VERIFY( (n == 0 || n == 3));

/* Computing MAX */
    d__4 = (d__1 = *a - *p, fabs(d__1)),	d__5 = (d__2 = *b - *f, fabs(d__2));
	d__4 = (d__4 >= d__5) ? d__4 : d__5,	d__5 = (d__3 = *h - *d, fabs(d__3));
    al = (d__4 >= d__5) ? d__4 : d__5;
/* Computing MIN */
    d__1 = *a + *p,	d__2 = *b + *f;
    d__1 = (d__1 <= d__2) ? d__1 : d__2,	d__2 = *h + *d;
    am = (d__1 <= d__2) ? d__1 : d__2;
    q = 0.;
    if (n == 3) q = .5;
    i__1 = (int) (am + q);
    for (k = (int) (al + q); k <= i__1; ++k) {
		ak = k - q;
		ret_val += (ak * 2. + 1.) *
					w_(a, p, d, h, &ak, g) *
					w_(b, f, h, d, &ak, e) *
					w_(a, p, b, f, &ak, c);
    }
    return ret_val;
} /* wnj_ */

/* ******************************************************************** */
/* FUNCTION THAT CALCULATES CLEBSCH-GORDAN COEFFICIENTS. THE FORMULAE   */
/* CAN BE FOUND ON PAGE 57 OF RICHARD N. ZARE'S `ANGULAR MOMENTUM'      */
/* ******************************************************************** */
double CG1(double j1, double j2, double j3,
			double p1, double p2, double p3)
{
    /* System generated locals */
    double ret_val, d__1;


    ret_val = 0.;
    if (j1 + j2 < j3 || (d__1 = j1 - j2, fabs(d__1)) > j3) return ret_val;
    if ( j1 < p1 || j2 < p2 || j3 < p3) return ret_val;
    if (p3 == p1 + 1. && p2 == 1.) {	 /* P' = P" + 1 */
        /* rR-BRANCH */
    	if (j3 == j1 + 1.) {
			ret_val = sqrt((j1 + p3) * (j1 + p3 + 1.) /
							(j1 * 2. + 1.) / (j1 * 2. + 2.));
    	}
		/* rQ-BRANCH */
    	else if (j3 == j1) {
			ret_val = -sqrt((j1 + p3) * (j1 - p3 + 1.) /
							(j1 * 2.) / (j1 + 1.));
    	}
		/* rP-BRANCH */
    	else if (j3 == j1 - 1.) {
			ret_val = sqrt((j1 - p3) * (j1 - p3 + 1.) /
							(j1 * 2.) / (j1 * 2. + 1.));
    	}
    }
    else if (p3 == p1 && p2 == 0.) {   /* P' = P" */
		/* qR-BRANCH */
    	if (j3 == j1 + 1.) {
			ret_val = sqrt((j1 - p3 + 1.) * (j1 + p3 + 1.) /
							(j1 * 2. + 1.) / (j1 + 1.));
    	}
		/* qQ-BRANCH */
    	else if (j3 == j1) {
			ret_val = p3 / sqrt(j1 * (j1 + 1.));
    	}
		/* qP-BRANCH */
    	else if (j3 == j1 - 1.) {
			ret_val = -sqrt((j1 - p3) * (j1 + p3) / j1 /
							(j1 * 2. + 1.));
    	}
    }
    else if (p3 == p1 - 1. && p2 == -1.) {	/* P' = P" - 1 */
		/* pR-BRANCH */
    	if (j3 == j1 + 1.) {
			ret_val = sqrt((j1 - p3) * (j1 - p3 + 1.) /
							(j1 * 2. + 1.) / (j1 * 2. + 2.));
    	}
		/* pQ-BRANCH */
    	else if (j3 == j1) {
			ret_val = sqrt((j1 - p3) * (j1 + p3 + 1.) /
							(j1 * 2.) / (j1 + 1.));
    	}
		/* pP-BRANCH */
    	else if (j3 == j1 - 1.) {
			ret_val = sqrt((j1 + p3 + 1.) * (j1 + p3) /
							(j1 * 2.) / (j1 * 2. + 1.));
    	}
    }
    return ret_val;
} /* cg1_ */

