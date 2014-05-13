#include "msdef.h"
#include "mathutil.h"
#include <malloc.h>
#include <math.h>
#include "model.h"
#include <complex>
using namespace std;
using std::complex;

/*** Local constant definitions ****/

#define N_CONST_LOCAL 11
//#define N_CONST_TWOFOLD 21
#define MBASE			1 /*** NO SPIN ***/
#define N_VOLATILE_TWF  4
#define N_EVSTATES      2
#define N_SYMM_STATES   3
#define M_NOISE         1e-6
#define NSORT			2
#define nConstGround    18
#define nConstExcited   30

#ifndef PI_CONST
#define PI_CONST 3.1415926526
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG 57.29577953
#endif

#define LOCK_CONSTANTS   1

typedef struct {
int N;
int K;
int WF;
} idx;

/*
  ***********************************************************************
  THIS MODEL IS USED FOR CALCULATION OF THE SPECTRA INVOLVING
  VIBRONIC FOURFOLD IN HUND'S CASE "B" BASIS SET IN THE BASIS SET OF
  LOCALIZED WAVEFUNCTIONS IN THREE EQUIVALENT POTENTIAL WELLS TRANSFORMED
  FROM AN A1, A2, E', E'' BASIS SET.
  INTERACTIONS (DO NOT) INCLUDE SPIN-ORBIT, CORIOLIS, OPTIONAL
  SPIN-ROTATIONAL INTERACTIONS AND CENTRIFUGAL DISTORTION PARAMETERS.

  TRANSITIONS BETWEEN LOWER STATE A2 LEVELS (MWC'S HAM) AND EXCITED
  STATE FOURFOLD STATES MAY BE CALCULATED.  ADDITIONALLY ALL
  TRANSITIONS THAT MAY BE CALCULATED BY MWC'S HAMILTONIAN MAY ALSO
  BE CALCULATED HERE.

  05-21-2013
  First debug tried with limited transitions calculable.---TC

  **********************************************************************
*/

  /***************************************************************************************************************************/
/***************************************************************************************************************************/
/* At this point, we define several "shortcut" functions that we will use to calculate intensities, basically, to un-clutter
the code and make the latter mode readable. We first recognize that in case "b" intensity calculations we don't need a
general expression for the 3J- and 6J-symbols, as the options are fairly limited. The values of these Wigner coefficients
are pre-calculated and written as three separate functions, corresponding to the product of 6J-symbol {J", J', 1},{N', N", 1/2}
3J-symbol, {N",1,N'}{K", q, -K'}, and the phase factor (-1)^{J+S-K'} for the three SPHERICAL TENSOR q components. These can
be readily combined to obtain x,y,z, or x, or+/-, z transition moments. NOTE: since J and S quantum numbers are not mixed,
they result in the common factor for all contributions to the intensity from a particular state and will be lost when TM is
squared, hence they are omitted. Only (-1)^K' matters.*/

double F(double J, double K)
{
	return sqrt((J - K)*(J + K + 1));
}

double TDM6J(int J2p, int Np, int J2pp, int Npp)
{
	double out;
	double Jp, Jpp;
	double phase, rNorm;

	Jp = ((double)J2p) / 2.0;
	Jpp = ((double)J2pp) / 2.0;

	rNorm = sqrt((Jp + 0.5)*(Jpp + 0.5)*(2 * Np + 1)*(2 * Npp + 1));

	/* check if thiangle relationships are fulfilled. Note: J2p and J2pp are J values, DOUBLED*/
	if (abs(Np - Npp) > 1) return 0.0;
	if (abs(J2p - J2pp) > 2) return 0.0;

	if (J2pp == 2 * Npp + 1 && J2p == 2 * Np + 1) /* case F"1 --> F'1 transition */
	{
		phase = ((Np - Npp) % 2) ? -1.0 : 1.0;
		out = 0.5*phase* sqrt((double)(Npp + Np)*(3 + Npp + Np));
	}
	else if (J2pp == 2 * Npp - 1 && J2p == 2 * Np - 1) /* case F"2 --> F'2 transition */
	{
		phase = ((Np - Npp) % 2) ? -1.0 : 1.0;
		out = (Np*Npp == 0) ? 0.0 :
			-0.5*phase*sqrt((double)(Np + Npp - 1)*(Np + Npp + 2));
	}
	else if (J2pp == 2 * Npp - 1 && J2p == 2 * Np + 1) /* case F"2 --> F'1 transition */
	{
		phase = ((Np - Npp) % 2) ? -1.0 : 1.0;
		out = (Npp == 0) ? 0 : -0.5 *phase *sqrt((double)(Np - Npp + 2)*(Npp - Np + 1));
	}
	else if (J2pp == 2 * Npp + 1 && J2p == 2 * Np - 1)	/* case F"1 --> F'2 transition */
	{
		phase = ((Np - Npp) % 2) ? -1.0 : 1.0;
		out = (Np == 0) ? 0 : -0.5*phase* sqrt((double)(Npp - Np + 2)*(Np - Npp + 1));
	}
	else out = 0.0;

	return (rNorm < M_NOISE) ? 0 : out / rNorm;
}

/* similarly, for the 3J-symbol*/
/* NOTE: this particular function calculates ONLY ONE SPECIFIC 3J-symbol, which
is relevant for the electronic dipole moment transition, i.e.
N"  1  N'
K"  q  -K'
This is not supposed to be an universal 3J-symbol calculator.
The purpose of this "plug" is to see if we can speed calculations up. If not, we'll
revert to mathutil.c functions. This function will return 0 if q is not -1, 0, or +1
*/
double TDM3J(int Np, int Kp, int Npp, int Kpp, int q)
{
	double out, rNorm, phase;

	if (Kpp + q - Kp != 0) return 0;
	if (abs(Np - Npp) > 1) return 0;
	if (abs(q)> 1) return 0;

	phase = ((Npp + Kpp) % 2) ? -1.0 : 1.0;

	if (Np == Npp) rNorm = sqrt((double)Npp*(1 + Npp)*(1 + 2 * Npp))*sqrt(2.0);
	else if (Np == Npp + 1) rNorm = sqrt((double)(1 + Npp)*(1 + 2 * Npp)*(1 + 2 * Np))*sqrt(2.0);
	else if (Np == Npp - 1) rNorm = sqrt((double)Npp*(1 + 2 * Np)*(1 + 2 * Npp))*sqrt(2.0);
	else rNorm = 0.0;

	if (Np == Npp)
		switch (q)
	{
		case 0: out = -phase*Kpp*sqrt(2.0);
			break;
		case 1: out = -phase*F(Npp, Kpp);
			break;
		case -1: out = phase*F(Npp, Kpp - 1);
			break;
		default: out = 0.0;
	}
	else if (Np == Npp + 1)
		switch (q)
	{
		case 0: out = -phase*sqrt((double)(Np*Np - Kpp*Kpp))*sqrt(2.0);
			break;
		case 1: out = phase*sqrt((double)(Np + Kpp)*(Np + Kpp + 1));
			break;
		case -1: out = phase*sqrt((double)(Np - Kpp)*(Np - Kpp + 1));
			break;
		default: out = 0.0;
	}
	else if (Np == Npp - 1)
		switch (q)
	{
		case 0: out = phase*sqrt((double)(Npp*Npp - Kpp*Kpp))*sqrt(2.0);
			break;
		case 1: out = phase*sqrt((double)(Np - Kpp)*(Np - Kpp + 1));
			break;
		case -1: out = phase*sqrt((double)(Np + Kpp)*(Np + Kpp + 1));
			break;
		default: out = 0.0;
	}
	return (rNorm < M_NOISE) ? 0.0 : out / rNorm;
}

/* functions MWC uses in Hamilt */
double G(double N, double S, double J)
{
	return N*(N + 1) + S*(S + 1) - J*(J + 1);
}

double C(double J, double N)
{
	double S = 0.5;
	return J*(J + 1) - N*(N + 1) - S*(S + 1);
}

double Fe(double J, double N)
{
	return -C(J, N) / 2 / N / (N + 1);
}

double Ro(double J, double N)
{
	double S = .5;
	return -(3 * Fe(J, N)*(C(J, N) + 1) + 2 * S*(S + 1)) / (2 * N - 1) / (2 * N + 3);
}

double P(double J, double N)
{
	double S = .5;
	return (N - J + S)*(N + J + S + 1);
}

double Q(double J, double N)
{
	double S = .5;
	return (J - N + S)*(N + J - S + 1);
}

double Psi(double J, double N)
{
	return -1. / N*sqrt(P(J, N)*Q(J, N - 1) / (2 * N - 1) / (2 * N + 1));
}

double g(double X, double Y)
{
	return sqrt((X - Y)*(X - Y - 1));
}

double JJNN(int Jp, int Jpp, int Np, int Npp)
{
	double result = sqrt((double)((2 * Jp + 1) * (2 * Jpp + 1) * (2 * Np + 1) * (2 * Npp + 1)));
	return result;
}
//MWC's functions for NSSW
double NSSW(int N, int K)
{
	int M;
	/*-------------Nuclear Spin Statistic Weight------------------*/
	M = K % 6;
	if (M == 1 || M == 5) /* K = 6n+1 */
	{
		return 0.0;
	}
	else if (M == 2 || M == 4) /* K = 6n+2
							   */
	{
		return 0.0;
	}
	else if (K == 0 && N % 2 == 0) /* KL = 0, NL = even */
	{
		return 0.0;
	}
	else
	{
		return 3.0;
	}

	/*-------------Nuclear Spin Statistic Weight------------------*/
}
UINT StatWeight(UINT nLoStateType, int* LoStQN, double* pdParam)
{
	int JL, NL, KL, SL, M, NSSW;
	JL = LoStQN[0];
	NL = LoStQN[1];
	KL = LoStQN[2];
	SL = (double)LoStQN[3] / 2.;
	int iLo = nLoStateType;
	if (nLoStateType == 0) /* iLo = 0 */
	{
		if (pdParam[9] == 1)  /* pdParam[9] = 1, NSSW is considered */
		{ /* MWC: nuclear spin statistical weight function for intensity calculation */
			M = ((int)(abs(KL))) % 6;
			if (M == 1 || M == -5) /* K = 6n+1 */
			{
				NSSW = 0;
			}
			else if (M == 2 || M == -4) /* K = 6n+2 */
			{
				NSSW = 0;
			}
			else if (KL == 0 && (int)(abs(NL)) % 2 == 0) /* KL = 0, NL = even */
			{
				NSSW = 0;
			}
			else
			{
				NSSW = 3 * (2 * JL + 1);
			}
		}
		if (pdParam[9] == 1)  /* pdParam[9] = 1, NSSW is not considered */
		{ /* MW: nuclear spin statistical weight function for intensity calculation */
			NSSW = 2 * JL + 1;
		}
		else
		{
			NSSW = 2 * JL + 1;
		}
	}
	return NSSW;
}

double lfabs(double x) { return (x<0.0) ? -x : x; }

/* Functions used in the Localized Hamiltonian (no spin) */
double Fpm(int N, int K)
{
	return (2 * N * (N + 1) - K * K);
}
double Fpp(int N, int K)
{
	return sqrt(N * (N + 1) - (K - 1) * (K - 2)) * sqrt(N * (N + 1) - K * (K - 1));
}
double Fmm(int N, int K)
{
	return sqrt(N * (N + 1) - (K + 1) * (K + 2)) * sqrt(N * (N + 1) - K * (K + 1));
}
double Fzz(int N, int K)
{
	return (K + 2) * (K + 1);
}
/* Model Name */
char* Name()
{
  return "NO3_Localized";
}

/*-------------------------State types-------------------------*/
/* Number of state types. The actual number of the state types is 2, however
by specifying A' and A" isolated states as separate types we can automatically
take their symmetry properties into account by selectively applying the components
of the transition moment in intensity calculations.*/
// HT: Working from here
UINT StaNum()
{
  return 1;
}

char* StName( UINT nStateType )
{
	switch (nStateType)
	{
	case 0:
		return "Localized";
	default:
		return "Localized";
	}
}

/*--------------------------Constants--------------------------*/
/* Number of constants in a given state type. A' and A" states have 10 constants each and
twofold has 2x10 plus magnitudes of SO and Cor parameters, tilt angle, total separation
Delta E = E(A')-E(A")*/
UINT ConNum( UINT nStateType )
{
	switch(nStateType)
	{
	case 0:
		return N_CONST_LOCAL;
	default : return 0;
	}

}

char* ConName( UINT nStateType, UINT nConst )
{
	switch(nStateType)
	{
	case 0:
		switch(nConst)
		{
			case 0 : return "C:A1";
			case 1 : return "B:A1";
			case 2 : return "C:A2";
			case 3 : return "B:A2";
			case 4 : return "C:E";
			case 5 : return "B:E";
			case 6 : return "h1:A1_E";
			case 7 : return "h1:A2_E";
			case 8 : return "h1:E_E";
			case 9 : return "dE1";
			case 10 : return "dE2";
			default : return "Constant";
		}//end switch nConst
	default : return "State"; //None of the defined
	}//end switch nStateType
}//end ConName

/* Initialize constants (defaults)*/
void InitCo( UINT nStateType, double* pdConst )
{
	for(int i = 0; i < (int)ConNum( nStateType ); i++)//initializes all constants to 0.0
	{
		pdConst[i] = 0.0;//I changed this from pdConst[i++]
	}

	switch (nStateType)
	{
	case 0 :
		pdConst[0] = 0.21;//C_A1
		pdConst[1] = 0.43;//B_A1
		pdConst[2] = 0.21;//C_A2
		pdConst[3] = 0.43;//B_A2
		pdConst[4] = 0.21;//C_E
		pdConst[5] = 0.43;//B_E
		return;
	}
}

/*--------------------------Parameters--------------------------*/
/* Number of parameters */
UINT ParNum()
{
  return 25;
}

/* Parameter names */
char* ParName( UINT nParam )//changed this to copy MWC
{
  switch( nParam ){
  case 7 : return "MaxJ";
  case 8 : return "MaxDltK";
  case 9 : return "parallell weighting";  //a-type
  case 10 : return "perpendicular weighting";  //b-type
  case 11 : return "null";  //c-type
  case 12 : return "LoSt K";
  case 13 : return "UpSt K";
  case 14 : return "LoSt K restrict";
  case 15 : return "UpSt K restrict";
  case 16 : return "LoSt N";
  case 17 : return "UpSt N";
  case 18 : return "LoSt N restrict";
  case 19 : return "UpSt N restrict";
  case 20 : return "P";
  case 21 : return "Q";
  case 22 : return "R";
  case 23 : return "NSSW";
  case 24 : return "Cmin";
  default : return "Parameter";/* Non of the defined*/
  }
}

/*Cmin is the minimum value of the expansion coefficient in the LOWER state for the purpose of
intensity calculations*/

/* Initialize parameters (defaults) */
void InitPa( double* pdParam )
{
  pdParam[0] = 10.0;    /* Parameters. T(K)*/
  pdParam[1] = 0.004;  /* Parameters. DelW(Dopp)*/
  pdParam[2] = 0.004;  /* Parameters. DelW(Norm)*/
  pdParam[3] = 0.001;    /* Parameters. Min Intensity*/
  pdParam[4] = 10.0;   /* Parameters. SpaceToSkip*/
  pdParam[5] = 0.0005; /* Parameters. TheorPlotRes*/
  pdParam[6] = 0.0;    /* Units 0.0<=>cm-1; 1.0<=>MHz; 2.0<=>GHz */
  pdParam[7] = 11;     /* Parameters. Jmax*/
  pdParam[8] = 5;      /* Maximum Delta K*/
  pdParam[9] = 1.0;     /* weight of a-type transition intensities.*/
  pdParam[10] = 0.0;   /* weight of b-typetransition intensities.*/
  pdParam[11] = .0;    /* weight of c-type transition intensities.*/
  pdParam[23] = 1.0;    /* NSSW */
  pdParam[24] = 0.0005;  /* whatever Dmitry's constant Cmin is */

  int i;
  for( i = 12; i < 20; pdParam[i++] = 0.0 );

  int j;
  for( j = 20; j < 23; pdParam[j++] = 1.0 );
}

/*------------------------Quantum Numbers------------------------*/
/* Number of "GOOD" Quantum Numbers -- in all models only J is conserved, so we leave it as is */
UINT QNnumG(UINT nStateType)
{
  return 1; /* J */
}

/* Number of "BAD" Quantum Numbers */
UINT QNnumB(UINT nStateType)
{
	switch(nStateType)
	{
		case 0 : return 3;//N, K, WF (N is good without spin, if this model does not include spin)
		default : return 3;
	}
}

/* Quantum Number names */
char* QNName( UINT nStateType, UINT nQNumb )
{
	switch(nStateType)
	{
	case 0 :
		switch(nQNumb)
		{
		case 0 : return "J";
		case 1 : return "N";
		case 2 : return "K";
		case 3 : return "WF";
		default : return "QN";/* Non of the defined*/
		}
	}
}

/* All Quantum Numbers are stored as integer, but we should know when
   they actually half integer. In that case we store doubled value and
   base is "2" */
UINT QNBase( UINT nStateType, UINT nQN ) //HT: Just returns 1
{
	switch(nStateType)
	{
	case 0 :
		switch (nQN)//j is half integer
		{
		case 0 : return MBASE;//Dmitry has this stored as 2
		case 3 : return MBASE;//for parity in MWC's model
		default : return 1;
		}
	}
}
/* Quantum number set that is first in the loop */
void StartQN( UINT nStateType, int* pnQNd )
{
  pnQNd[0]=0; /* HT: This corresponds to J=N=0 */
}

/* The quantum number loop iteration (continue) condition. The Jmax value is specified in
as integer, which is no big deal since it only sets the limit, but for comparison purpose
we need to convert it to doubled form.*/
BOOL ContQN( UINT nStateType, int* pnQNd, double* pdParam )
{
  return(pnQNd[0] <= 2*pdParam[7]);//changed from < to <= to match MWC, I think he's right
}
//HT: Is this <=pdParam[7] with MBASE=1?

/* The quantum number loop iteration(next set in the loop). */
void NextQN( UINT nStateType, int* pnQNd, double* pdParam )
{
  pnQNd[0]+=1;
  return;
}

/* OK. As you could have noticed that the twofold model is fairly involved in terms of the quantum number
indexing, ESPECIALLY in case "b" where J+1/2 and J-1/2 subblocks are not of the same size. Therefore,
to facilitate calculation of the matrix elements and inserting them into the matrix properly, we will
write two book-keeping functions, one of which translates the linear array index i or j into the
set of quantum numbers that the row or column corresponds to, and the other function does the opposite
thing, i.e. identifies i or j by the set of quantum numbers. */

idx lIndex(int *pnQNd, int nx)
{
	int Ks, spot;
	idx out;

	Ks = 2*pnQNd[0]+1; //Number of K values in each WF block
	int sBlock = nx / Ks;
	out.WF = sBlock;

	spot = (nx + Ks) % Ks;
	out.K = spot - (Ks - 1) / 2;

	out.N = pnQNd[0];

	return out;
}

/*The second indexing function is for the reverse indexing.
Here, we don't check for the validity of quantum numbers, i.e.
nonnegative N, etc.*/

int lReverseIndex(int *pnQNd, int N, int K, int WF)
{
  int Ks = (2 * pnQNd[0]) + 1;
  int nx = Ks * WF + K + pnQNd[0];
  return nx;
}

/*------------------------Main Part------------------------*/
/* Assign() function is not required for the model operation but it helps to give meaningful
assignments. IMPORTANT NOTE: in this model assigment is given in terms of Wang components,
therefore K assumes only non-negative values. If such an assignment is implemented, DO NOT
use the MaxDltK option at all, since this will confuse the core and will result in
missing transitions. Also, in this program the symmetry labels will be given to the
eigenvectors, which currently are not used anywhere except in the user-end output.*/

/* NOTE 2: one of the issues is that the mechanism implemented in Specview deals with one
eigenvector at a time, whereas the algorithm implemented in Mathematica version and in
my earlier program for fitting methoxy analyzes the entire J-subblock preventing
duplicate assignments. In weakly coupled systems this won't be the issue, however it is
not clear how the assignment (labeling) function will behave in near-degenerate cases. */

void Assign( UINT nStateType, double* pdStWF, int* StQN, int StNum )
{
	if(nStateType == 0)
	{
		idx QN;
		int id, DIM;
		DIM = 4*(StQN[0] + 1);
		double max = abs(pdStWF[0]);
		id = 0;
		for(int i = 1; i < DIM; i++)
		{
			if(abs(pdStWF[i]) > max)
			{
				max = abs(pdStWF[i]);
				id = i;
			}
		}
		QN = lIndex(StQN, id);
		StQN[1] = QN.K;
		StQN[2] = QN.WF;
		return;
	}//My simple assign function
}//end assign

/* Get Hamiltonian matrix size for a given J block */
UINT HamSize( UINT nStateType, int* pnQNd )
{
  return 4*(2*pnQNd[0] + 1);  // HT: Symm * (2*J+1)
}


/* Building Hamiltonian blocks for a given J, includes code for all supported states*/
void Hamilt( UINT nStateType,
	     double* pdConst,
	     int* pnQNd,
	     std::complex<double>** ppdH,int** ppnQNm)
{
	if(nStateType == 0) 
	{
	  int N,K,WF,Ks,DIM;
	  idx QN;
	  double J;
	  double CA1, BA1, CA2, BA2, CE, BE, h1A1E, h1A2E, h1EE, E1, E2;
	  std::complex<double> G1o2,Gn1o6,G7o6,Gn2o3,G2o3,Gn1o3,G1o3,Gn1o1;
	  std::complex<double> GC1o2,GCn1o6,GC7o6,GCn2o3,GC2o3,GCn1o3,GC1o3,GCn1o1;
	  CA1 = pdConst[0]; 
	  BA1 = pdConst[1]; 
	  CA2 = pdConst[2]; 
	  BA2 = pdConst[3]; 
	  CE = pdConst[4]; 
	  BE = pdConst[5]; 
	  h1A1E = pdConst[6]; 
	  h1A2E = pdConst[7]; 
	  h1EE = pdConst[8]; 
	  E1 = pdConst[9]; 
	  E2 = pdConst[10]; 

	  /*Ghost Operator*/
	  /* G(host)(C)(onjugate)(n)(egative)Xo(ver)Y */
	  /* Parts in N-^2 */
	  std::complex<double> A1o2(0, PI_CONST/2);
	  G1o2=exp(A1o2);
	  std::complex<double> An1o6(0, -PI_CONST/6);
	  Gn1o6=exp(An1o6);
	  std::complex<double> A7o6(0, 7*PI_CONST/6);
	  G7o6=exp(A7o6);
	  std::complex<double> An2o3(0, -2*PI_CONST/3);
	  Gn2o3=exp(An2o3);
	  std::complex<double> A2o3(0, 2*PI_CONST/3);
	  G2o3=exp(A2o3);
	  std::complex<double> An1o3(0, -1*PI_CONST/3);
	  Gn1o3=exp(An1o3);
	  std::complex<double> A1o3(0, PI_CONST/3);
	  G1o3=exp(A1o3);
	  std::complex<double> An1o1(0, -1*PI_CONST);
	  Gn1o1=exp(An1o1);

	  /*Parts in N+^2*/
	  std::complex<double> AC1o2(0, -PI_CONST/2);
	  GC1o2=exp(AC1o2);
	  std::complex<double> ACn1o6(0, PI_CONST/6);
	  GCn1o6=exp(ACn1o6);
	  std::complex<double> AC7o6(0, -7*PI_CONST/6);
	  GC7o6=exp(AC7o6);
	  std::complex<double> ACn2o3 = A2o3;
	  GCn2o3=exp(ACn2o3);
	  std::complex<double> AC2o3 = An2o3;
	  GC2o3=exp(AC2o3);
	  std::complex<double> ACn1o3 = A1o3;
	  GCn1o3=exp(ACn1o3);
	  std::complex<double> AC1o3 = An1o3;
	  GC1o3=exp(AC1o3);
	  std::complex<double> ACn1o1(0, PI_CONST);
	  GCn1o1=exp(ACn1o1);

		J = pnQNd[0];
		Ks = 2 * J + 1; //Size of blocks
		DIM = 4 * Ks;

		for(int ii = 0; ii < DIM; ii++)  /* INITIALIZE THE HAMILTONIAN MATRIX*/
		{
			for(int jjj = 0; jjj < DIM; jjj++)
			{
				ppdH[ii][jjj] = 0.0;
			}
		}
		for(int ii=0; ii < DIM; ii++) /* ELEMENTS INDEPENDENT OF N */
		{
			//Get Quantum Numbers
			QN = lIndex(pnQNd, ii);

			K = QN.K;
			WF = QN.WF;
			N = QN.N;
			ppnQNm[ii][0] = N;
			ppnQNm[ii][1] = K;
			ppnQNm[ii][2] = WF;

			if(WF==0)
			{
				for(int jjj=0; jjj < DIM; jjj++)
				{
					//Get Quantum Numbers
					QN = lIndex(pnQNd, jjj);

					K = QN.K;
					WF = QN.WF;
					N = QN.N;
					ppnQNm[jjj][0] = N;
					ppnQNm[jjj][1] = K;
					ppnQNm[jjj][2] = WF;

					if(WF==0)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==2)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==3)
					{
						ppdH[ii][jjj] = E1/3;
					}
				}//End jjj loop
			}//End WF = 0
			if(WF==2)
			{
				for(int jjj=0; jjj < DIM; jjj++)
				{
					//Get Quantum Numbers
					QN = lIndex(pnQNd, jjj);

					K = QN.K;
					WF = QN.WF;
					N = QN.N;
					ppnQNm[jjj][0] = N;
					ppnQNm[jjj][1] = K;
					ppnQNm[jjj][2] = WF;

					if(WF==0)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==2)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==3)
					{
						ppdH[ii][jjj] = E1/3;
					}
				}//End jjj loop
			}//End WF = 2
			if(WF==3)
			{
				for(int jjj=0; jjj < DIM; jjj++)
				{
					//Get Quantum Numbers
					QN = lIndex(pnQNd, jjj);

					K = QN.K;
					WF = QN.WF;
					N = QN.N;
					ppnQNm[jjj][0] = N;
					ppnQNm[jjj][1] = K;
					ppnQNm[jjj][2] = WF;

					if(WF==0)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==2)
					{
						ppdH[ii][jjj] = E1/3;
					}
					if(WF==3)
					{
						ppdH[ii][jjj] = E1/3;
					}
				}//End jjj loop
			}//End WF = 3
		}//End ii Loop

		for(int ii=0; ii < DIM; ii++) /* ELEMENTS INDEPENDENT OF N (E2) */
		{
			//Get Quantum Numbers
			QN = lIndex(pnQNd, ii);

			K = QN.K;
			WF = QN.WF;
			N = QN.N;
			ppnQNm[ii][0] = N;
			ppnQNm[ii][1] = K;
			ppnQNm[ii][2] = WF;

			if(WF==1)
			{
				for(int jjj=Ks; jjj < 2*Ks; jjj++)
				{
					ppdH[ii][jjj] = E2;
				}
			}//End WF = 1
		}//End ii loop

		for(int i = 0; i < DIM; i++) /*START BUILDING HAMILTONIAN*/
		{
			//First, get Quantum Numbers
			//consider inlining this function call
			QN = lIndex(pnQNd, i);

			K = QN.K;
			WF = QN.WF;
			N = QN.N;
			ppnQNm[i][0] = N;
			ppnQNm[i][1] = K;
			ppnQNm[i][2] = WF;

/* Terms Linear in [N+,N-] */

			if(WF==0)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K, o);
					if(o==0) // S1
					{
						ppdH[i][j] += ((BA1 + 2 * BE) / 3) * Fpm(N, K);
					}
					if(o==2) // S2
					{
						ppdH[i][j] += ((BA1 - BE) / 3) * Fpm(N, K);
						ppdH[j][i] = ppdH[i][j];
					}
					if(o==3) // S3	
					{
						ppdH[i][j] += ((BA1 - BE) / 3) * Fpm(N, K);
						ppdH[j][i] = ppdH[i][j];
					}
				}// End o loop
			} //End WF==0
			if(WF==1)
			{
				ppdH[i][i] += BA2 * (2 * N * (N + 1) - K * K);
			} //End WF==1
			if(WF==2)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K, o);
					if(o==2) // S2
					{
						ppdH[i][j] += ((BA1 + 2 * BE) / 3) * Fpm(N, K);
					}
					if(o==3) // S3	
					{
						ppdH[i][j] += ((BA1 - BE) / 3) * Fpm(N, K);
						ppdH[j][i] = ppdH[i][j];
					}
				}// End o loop
			} //End WF==2
			if(WF==3)
			{
				ppdH[i][i] += ((BA1 + 2 * BE) / 3) * Fpm(N, K);

			} //End WF==3

/* Terms linear in N+^2 */

			if(WF==0) // S1
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K-2, o);
					if(j >= o * Ks) // Element is in intended WF block
					{
						if(o==0) // S1
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==1) // A2
						{
							ppdH[i][j] += h1A2E / sqrt(2) * GC1o2 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==2) // S2
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * GCn1o3 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * GC1o3 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==0
			if(WF==1) // A2
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K-2, o);
					if(j >= o * Ks) // Element is in intended WF block
					{
						if(o==2) // S2
						{
							ppdH[i][j] += (h1A2E / sqrt(2)) * GCn1o6 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += (h1A2E / sqrt(2)) * GC7o6 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==1
			if(WF==2) // S2
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K-2, o);
					if(j >= o * Ks) // Element is in intended WF block
					{
						if(o==2) // S2
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * GCn2o3 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * GCn1o1 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==2
			if(WF==3) // S3
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K-2, o);
					if(j >= o * Ks) // Element is in intended WF block
					{
						if(o==3) // S3
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * GC2o3 * Fpp(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==3

/* Terms Linear in N-^2 and Nz^2 (Respectively) */

			if(WF==0)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K+2, o);
					if(j < (o + 1) * Ks) // Element is in intended WF block
					{
						if(o==0) // S1
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * Fmm(N, K)  // N-^2
								+ (CA1 + 2 * CE) / 3 * Fzz(N, K); // Nz^2
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==1) // A2
						{
							ppdH[i][j] += h1A2E / sqrt(2) * G1o2 * Fmm(N, K)
								;
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==2) // S2
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * Gn1o3 * Fmm(N, K)
								+ (CA1 - CE) * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * G1o3 * Fmm(N, K)
								+ (CA1 - CE) * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==0
			if(WF==1)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K+2, o);
					if(j < (o + 1) * Ks) // Element is in intended WF block
					{
						if (o == 1) // A2
						{
							ppdH[i][j] +=
								CA2 * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==2) // S2
						{
							ppdH[i][j] += (h1A2E / sqrt(2)) * Gn1o6 * Fmm(N, K)
								;
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += (h1A2E / sqrt(2)) * G7o6 * Fmm(N, K)
								;
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==1
			if(WF==2)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K+2, o);
					if(j < (o + 1) * Ks) // Element is in intended WF block
					{
						if(o==2) // S2
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * Gn2o3 * Fmm(N, K)
								+ (CA1 + 2 * CE) / 3 * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
						if(o==3) // S3
						{
							ppdH[i][j] += ((h1A1E - h1EE) / 3) * Gn1o1 * Fmm(N, K)
								+ (CA1 - CE) * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==2
			if(WF==3)
			{
				for(int o = WF; o < 4; o++)
				{
					int j=lReverseIndex(pnQNd, N, K+2, o);
					if(j < (o + 1) * Ks) // Element is in intended WF block
					{
						if(o==3) // S3
						{
							ppdH[i][j] += (2 * h1A1E + h1EE) / 3 * G2o3 * Fmm(N, K)
								+ (CA1 + 2 * CE) / 3 * Fzz(N, K);
							ppdH[j][i] = ppdH[i][j];
						}
					} // End j conditional
				}// End o loop
			} //End WF==3
		}//end Hamiltonian for loop
	}//end if StateType == 0
	return;
}//end Hamilt



/* Get derivative of the Hamiltonian matrix on constant number k for a
   given good q.num. set */
void DHamilt( UINT nStateType,
              double* pdConst,
              int* pnQNd,
              std::complex<double> ** ppdDH, int k )
{

	  UINT i, j;
	  double *pdIndicators;
	  std::complex<double> **ppdHbase;
	  int **ppnQNm;
	  UINT nDim;
	  int lCN;
	  lCN = ConNum(nStateType);
	  pdIndicators = (double *)calloc(lCN , sizeof( double ) );

	  nDim = HamSize(nStateType,pnQNd); 
	  for ( int i = 0; i < lCN; i++ ) pdIndicators[i]=1.0;
	  pdIndicators[k] +=1.0;
	  ppnQNm = (int**)calloc( nDim, sizeof( int* ) );
	  ppdHbase = (std::complex<double>**)calloc(nDim, sizeof(double*));

	  for ( i = 0; i < nDim; i++ )
	  {
		ppnQNm[i] = (int*)calloc( QNnumB(nStateType), sizeof(int) );
		ppdHbase[i] = (std::complex<double>*)calloc( nDim, sizeof(double));
	  }
	Hamilt (nStateType, pdConst, pnQNd, ppdHbase, ppnQNm);
	Hamilt( nStateType, pdIndicators, pnQNd, ppdDH, ppnQNm);

	for (i=0; i< nDim; i++)
		for(j=0; j< nDim; j++)
			ppdDH[i][j] -= ppdHbase[i][j];

	for ( i = 0; i < nDim; i++ )
	{
	 free (ppdHbase[i]);
	 free( ppnQNm[i] );
	}
	free(ppnQNm);
	free (ppdHbase);
	free( pdIndicators );
}

//C.37, Dmitry's write up
double perpIntensity(int i, int* UpStQN, int* LoStQN, idx LoStIndex, int symm, int UpStateK, int UpStateN, double* pdLoStWF, double* pdUpStWF, double dWeightB)
{
      int j = lReverseIndex(UpStQN, UpStateN, UpStateK, symm);
      double coeff = pow(-1.0, (double)(UpStQN[0] + 1) / 2.0 + (double)UpStateK) / sqrt(2.0);
      double SixJ = TDM6J(UpStQN[0], UpStateN, LoStQN[0], LoStIndex.N);
      double ThreeJ_One = TDM3J(UpStateN, -1 * UpStateK, LoStIndex.N, LoStIndex.K, 1);
      double ThreeJ_Two = TDM3J(UpStateN, UpStateK, LoStIndex.N, LoStIndex.K, -1);
      double jjnn = JJNN(UpStQN[0], LoStQN[0], UpStateN, LoStIndex.N);
      double coeff_Two = pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1.0);//coefficient in front of second 3J symbol, p = symm - 1, k = j = 2.0
	  //Note to Terrance: Shouldn't these be upperstate constants?
      double co = pow(-1.0, LoStIndex.N + UpStateN + 1);//takes care of 3J symbol being different from function
      double intensity = jjnn * coeff * SixJ * co * (ThreeJ_One + coeff_Two * ThreeJ_Two) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
      return intensity;
  }//end perpIntensity function

//C.32, Dmitry's write up
double parIntensity(int i, int* UpStQN, int* LoStQN, idx LoStIndex, int UpStateN, double* pdLoStWF, double* pdUpStWF, double dWeightA)
{
    int j = lReverseIndex(UpStQN, UpStateN, LoStIndex.K, 0);//zero for symm because all are A1 levels
    double coeff = -1.0 * pow(-1.0, (double)LoStIndex.K + ((double)UpStQN[0] + 1.0) / 2.0);
    double co = pow(-1.0, LoStIndex.N + UpStateN + 1);
    double jjnn = JJNN(UpStQN[0], LoStQN[0], UpStateN, LoStIndex.N);
    double SixJ = TDM6J(UpStQN[0], UpStateN, LoStQN[0], LoStIndex.N);
    double ThreeJ = TDM3J(UpStateN, LoStIndex.K, LoStIndex.N, LoStIndex.K, 0);
    double intensity = coeff * jjnn * SixJ * ThreeJ * pdLoStWF[i] * pdUpStWF[j] * dWeightA;
    return intensity;
}//end parIntensity function

/* Intensity of a given transition */
void Inten( int* pnCount,
	    BOOL* pbCont,
	    double* pdParam,
	    double dLoStEnerg,
	    UINT nLoStateType,
	    double* pdLoStWF, int* LoStQN,
		double dUpStEnerg,
		UINT nUpStateType,
	    double* pdUpStWF, int* UpStQN,
	    float* fInten )
{

    if(nLoStateType == 0 && nUpStateType == 0)
    {
		double dNorm, sum, BF, LS;
		int f, dJL, dJU, DIM, uNL, uNH, j, coeff, symm;
		double static dWeightA, dWeightB, dWeightCsq, EMIN, RK, RINTE;
		BOOL static bFlagA, bFlagBC;
		//From model.c
		  if (!*pbCont)
		  {/* This Ground Level was the first one, so you should*/
			  /* setup trap for various ladders. Keep in mind that */
			  EMIN = dLoStEnerg; /* states are sorted and lowest*/
			  *pbCont = TRUE;     /* are coming first. Here we have just one*/
			  /* BOLTZMANN FACTOR hv/kt-->v/(k/h)t, k/h in Hz/K */
			  RK = 1.380662E-23 / 6.626176E-34;
			  RK = (pdParam[6] == 2.0)?RK * 1E-9:     /* GHz*/
				  ((pdParam[6] == 1.0)?RK * 1E-6:RK * 1E-9 / 29.9792458 );
					/*MHz               cm-1*/
			  dNorm = fabs(pdParam[9]) + fabs(pdParam[10]) + fabs(pdParam[11]);
			  dNorm = (dNorm==0.0)?1.0:dNorm;
			  dWeightA = sqrt(fabs(pdParam[9])/dNorm);     /* weight of a-type transition intensities.*/
			  dWeightB = sqrt(fabs(pdParam[10]/dNorm));   /* weight of b-typetransition intensities.*/
			  dWeightCsq = fabs(pdParam[11]/dNorm);    /* weight of c-type transition intensities.*/
			  bFlagA =  (dWeightA != 0.0);        /*ladder.*/
			  bFlagBC = ( dWeightB!=0.0 || dWeightCsq!=0.0);
		  }
		  *fInten = 0.0f;
		  sum = 0.0;
		  dJL = *LoStQN;
		  dJU = *UpStQN;
		  f = abs(dJL - dJU);
		  if(f > 3)//checks delta J = 0,+/- 1
		  {
			  return;
		  }
		  //DIM = 2*(LoStQN[0] + 1);//sets the bounds for the for loop over the GS EV components
		  DIM = HamSize(1, LoStQN);
		  idx LoStIndex;
		  //First find out what N values are in the upper state based on it's J to make sure we don't get erroneous indices from lReverseIndex
		  uNL = (UpStQN[0]-1)/2;
		  uNH = (UpStQN[0]+1)/2;
		  int upStateN;
		  int upStateK;

		  if(bFlagA)//C.32
		  {
			  for(int i = 0; i < DIM; i++)
			  {
				  LoStIndex = lIndex(LoStQN, i);// MWCIndex(LoStQN, i);
				 //Now check for 3 possible components that the GS component can link with first checking to see if the bad QN's exist in the excited state and only check for A1 (Symm = 0)
				  if(LoStIndex.N - 1 == uNL || LoStIndex.N - 1 == uNH)
				  {
				      upStateN = LoStIndex.N - 1;
				      sum += parIntensity(i, UpStQN, LoStQN, LoStIndex, upStateN, pdLoStWF, pdUpStWF, dWeightA);
                  }
				  if((LoStIndex.N == uNL || LoStIndex.N == uNH) && LoStIndex.K != 0 && LoStIndex.N != 0)
				  {
				      upStateN = LoStIndex.N;
				      sum += parIntensity(i, UpStQN, LoStQN, LoStIndex, upStateN, pdLoStWF, pdUpStWF, dWeightA);
                  }
				  if(LoStIndex.N + 1 == uNL || LoStIndex.N + 1 == uNH)
				  {
				      upStateN = LoStIndex.N + 1;
				      sum += parIntensity(i, UpStQN, LoStQN, LoStIndex, upStateN, pdLoStWF, pdUpStWF, dWeightA);
                  }
			  }//end for loop
		  }//end if bFlagA == true

		  if(bFlagBC)//C.37
		  {
			//  for(int symm = 2; symm < 4; symm++)
			//  {
				  for(int i = 0; i < DIM; i++)
				  {
				      LoStIndex = lIndex(LoStQN, i);//MWCIndex(LoStQN, i);
				      if(LoStIndex.N - 1 == uNL || LoStIndex.N - 1 == uNH)
					  {
					      upStateN = LoStIndex.N - 1;
						  if(LoStIndex.K - 1 >= upStateN * -1)
						  {
						      upStateK = LoStIndex.K - 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
						  if(LoStIndex.K + 1 <= upStateN)//I think this is the key!!!!!!!!!!!!
						  {
						      upStateK = LoStIndex.K + 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
					  }//end delta N = -1

					  if((LoStIndex.N == uNL || LoStIndex.N == uNH) && LoStIndex.N != 0)
					  {
					      upStateN = LoStIndex.N;
						  if(LoStIndex.K - 1 >= upStateN * -1)
						  {
						      upStateK = LoStIndex.K - 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
						  if(LoStIndex.K + 1 <= upStateN)
						  {
						      upStateK = LoStIndex.N + 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
					  }//end delta N = 0

					  if(LoStIndex.N + 1 == uNL || LoStIndex.N + 1 == uNH)
					  {
					      upStateN = LoStIndex.N + 1;
						  if(LoStIndex.K - 1 >= upStateN * -1)
						  {
						      upStateK = LoStIndex.K - 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
						  if(LoStIndex.K + 1 <= upStateN)
						  {
						      upStateK = LoStIndex.K + 1;
						      sum += perpIntensity(i, UpStQN, LoStQN, LoStIndex, symm, upStateK, upStateN, pdLoStWF, pdUpStWF, dWeightB);
						  }
					  }//end delta N = +1
				  }//end loop over lostate wf
			 // }//end loop over E+/E-
		  }//end if bFlagBC

		  LS = sum*sum; /* "Line strength" multiplies square of the MAG(nitude) by (2J'+1)(2J"+1) */
		  BF = exp(-(dLoStEnerg - EMIN)/(RK*pdParam[0]))-exp(-(dUpStEnerg - EMIN)/(RK*pdParam[0]));   /* Boltzmann factor */
		  RINTE = ((int)pdParam[23]==0? 1 : NSSW(LoStQN[1], LoStQN[2])) * BF * LS;
		  *(pnCount)++;
		  if((RINTE >= pdParam[3]))// && (dUpStEnerg>dLoStEnerg))
		  {
			  *fInten = (float)RINTE;
		  }
		  return;
    }//end upstate == lowstate == 0

  }//end Inten function
  //}
