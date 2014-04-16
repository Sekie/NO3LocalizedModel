#include "msdef.h"
#include "mathutil.h"
#include <malloc.h>
#include <math.h>
#include "model.h"
using namespace std;


/*** Local constant definitions ****/

#define N_CONST_SINGLE	18
#define N_CONST_TWOFOLD 26
#define MBASE			2
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
int Symm;
int N;
int K;
int P;
} idx;

/*
  ***********************************************************************
  THIS MODEL IS USED FOR CALCULATION OF THE SPECTRA INVOLVING
  VIBRONIC FOURFOLD IN HUND'S CASE "B" BASIS SET AND A1, A2, E+ and E-
  VIBRONIC BASIS SET. INTERACTIONS INCLUDE SPIN-ORBIT, CORIOLIS, OPTIONAL
  SPIN-ROTATIONAL INTERACTIONS AND CENTRIFUGAL DISTORTION PARAMETERS.

  THIS MODEL INCLUDES MING-WEI CHEN'S (MWC) OBLATE SYMMETRIC TOP HAMILTONIAN
  WITH SPIN ROTATION AND CENTRIFUGAL DISTORTION AS A SUBSET. HIS ORIGINAL CODE
  IS USED WITHOUT MODIFICATION WITH HIS PERMISSION. WHERE HIS CODE IS USED IT IS
  LABELED AS 'MWC'.

  TRANSITIONS BETWEEN LOWER STATE A2 LEVELS (MWC'S HAM) AND EXCITED
  STATE FOURFOLD STATES MAY BE CALCULATED.  ADDITIONALLY ALL
  TRANSITIONS THAT MAY BE CALCULATED BY MWC'S HAMILTONIAN MAY ALSO
  BE CALCULATED HERE.

  MANY OF THE SHORTCUT FUNCTIONS AND EXPLANATORY COMMENTS WERE WRITTEN BY
  DMITRY MELNIK. THESE ARE LEFT IN PLACE BECAUSE IT WAS HIS TWO-FOLD MODEL
  THAT WAS USED AS A TEMPLATE FOR THE FOURFOLD MODEL.
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


double F( double J, double K )
{
  return sqrt((J-K)*(J+K+1));
}

double TDM6J(int J2p, int Np, int J2pp, int Npp)
{
double out;
double Jp, Jpp;
double phase, rNorm;

Jp = ((double)J2p)/2.0;
Jpp = ((double)J2pp)/2.0;

rNorm = sqrt((Jp+0.5)*(Jpp+0.5)*(2*Np+1)*(2*Npp+1));

/* check if thiangle relationships are fulfilled. Note: J2p and J2pp are J values, DOUBLED*/
if (abs(Np-Npp) > 1) return 0.0;
if (abs(J2p-J2pp) > 2 ) return 0.0;

if ( J2pp == 2*Npp+1 && J2p == 2*Np+1 ) /* case F"1 --> F'1 transition */
{
	phase = ((Np -Npp)%2)? -1.0:1.0;
	out = 0.5*phase* sqrt((double)(Npp+Np)*(3+Npp+Np));
}
else if ( J2pp == 2*Npp-1 && J2p == 2*Np-1 ) /* case F"2 --> F'2 transition */
{
	phase = ((Np -Npp)%2)? -1.0:1.0;
	out = (Np*Npp==0)? 0.0:
		-0.5*phase*sqrt((double)(Np+Npp-1)*(Np+Npp+2));
}
else if ( J2pp == 2*Npp-1 && J2p == 2*Np +1) /* case F"2 --> F'1 transition */
{
	phase = ((Np -Npp)%2)? -1.0:1.0;
	out = (Npp==0)? 0: -0.5 *phase *sqrt((double)(Np-Npp+2)*(Npp-Np+1));
}
else if (J2pp == 2*Npp+1 && J2p == 2*Np -1)	/* case F"1 --> F'2 transition */
{
	phase = ((Np -Npp)%2)? -1.0:1.0;
	out = (Np==0)? 0: -0.5*phase* sqrt((double)(Npp-Np+2)*(Np-Npp+1));
}
else out = 0.0;

return (rNorm < M_NOISE )? 0: out/rNorm;
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

if (Kpp + q -Kp !=0) return 0;
if (abs(Np - Npp ) > 1) return 0;
if (abs(q)> 1) return 0;

phase = ((Npp+Kpp)%2)? -1.0: 1.0;

if (Np==Npp) rNorm = sqrt((double)Npp*(1+Npp)*(1+2*Npp))*sqrt(2.0);
else if (Np == Npp+1) rNorm = sqrt((double)(1+Npp)*(1+2*Npp)*(1+2*Np))*sqrt(2.0);
else if (Np == Npp-1) rNorm = sqrt((double)Npp*(1+2*Np)*(1+2*Npp))*sqrt(2.0);
else rNorm =0.0;

if (Np==Npp)
	switch (q)
	{
	case 0: out = -phase*Kpp*sqrt(2.0);
		break;
	case 1: out = -phase*F(Npp, Kpp);
		break;
	case -1: out = phase*F(Npp,Kpp-1);
		break;
	default: out = 0.0;
	}
else if (Np==Npp+1)
	switch (q)
   {
	case 0: out = -phase*sqrt((double)(Np*Np - Kpp*Kpp))*sqrt(2.0);
		break;
	case 1: out = phase*sqrt((double)(Np+Kpp)*(Np+Kpp+1));
		break;
   case -1: out = phase*sqrt((double)(Np-Kpp)*(Np-Kpp+1));
		break;
	default: out = 0.0;
	}
else if (Np==Npp-1)
	switch (q)
	{
	case 0: out = phase*sqrt((double)(Npp*Npp - Kpp*Kpp))*sqrt(2.0);
		break;
	case 1: out = phase*sqrt((double)(Np-Kpp)*(Np-Kpp+1));
		break;
	case -1: out = phase*sqrt((double)(Np+Kpp)*(Np+Kpp+1));
		break;
	default: out = 0.0;
	}
return (rNorm < M_NOISE)? 0.0: out/rNorm;
}

/* functions MWC uses in Hamilt */
double G(double N, double S, double J)
{
	return N*(N+1) + S*(S+1) - J*(J+1);
}

double C(double J, double N)
{
	double S = 0.5;
	return J*(J+1) - N*(N+1) - S*(S+1);
}

double Fe(double J, double N)
{
	return -C(J, N)/2/N/(N+1);
}

double Ro(double J, double N)
{
	double S = .5;
	return -(3*Fe(J, N)*(C(J,N)+1)+2*S*(S+1))/(2*N-1)/(2*N+3);
}

double P(double J, double N)
{
	double S = .5;
	return (N-J+S)*(N+J+S+1);
}

double Q(double J, double N)
{
	double S = .5;
	return (J-N+S)*(N+J-S+1);
}

double Psi(double J, double N)
{
	return -1./N*sqrt(P(J,N)*Q(J,N-1)/(2*N-1)/(2*N+1));
}

double g(double X, double Y)
{
	return sqrt((X-Y)*(X-Y-1));
}

double JJNN(int Jp, int Jpp, int Np, int Npp)
{
	double result = sqrt((double)((2 * Jp + 1) * (2 * Jpp + 1) * (2 * Np + 1) * (2 * Npp + 1)));
	return result;
}
//MWC's functions for NSSW
double NSSW(int N, int K)
{ int M ;
    /*-------------Nuclear Spin Statistic Weight------------------*/
 				M = K%6;
 				if ( M == 1 || M == 5 ) /* K = 6n+1 */
 				{
 					return 0.0;
 				}
 				else if ( M == 2 || M == 4 ) /* K = 6n+2
*/
 				{
 					return 0.0;
 				}
 				else if ( K == 0 && N %2 == 0 ) /* KL = 0, NL = even */
 				{
 					return 0.0;
 				}
 				else
 				{
 					return 3.0 ;
 				}

    /*-------------Nuclear Spin Statistic Weight------------------*/
}
UINT StatWeight(UINT nLoStateType, int* LoStQN, double* pdParam)
{
  int JL, NL, KL, SL, M, NSSW;
  JL= LoStQN[0];
  NL= LoStQN[1];
  KL= LoStQN[2];
  SL=(double)LoStQN[3]/2.;
  int iLo = nLoStateType;
  if (nLoStateType == 0) /* iLo = 0 */
	{
		if (pdParam[9] == 1)  /* pdParam[9] = 1, NSSW is considered */
		{ /* MWC: nuclear spin statistical weight function for intensity calculation */
				M = ((int)(abs(KL)))%6;
				if ( M == 1 || M == -5 ) /* K = 6n+1 */
				{
					NSSW = 0;
				}
				else if ( M == 2 || M == -4 ) /* K = 6n+2 */
				{
					NSSW = 0;
				}
				else if ( KL == 0 && (int)(abs(NL))%2 == 0 ) /* KL = 0, NL = even */
				{
					NSSW = 0;
				}
				else
				{
					NSSW = 3 * ( 2 * JL + 1 );
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

double lfabs(double x) {return (x<0.0)? -x:x;}

//functions I use in my Hamiltonian
double FNK(int N, int K)
{
	double n = N;
	double k = K;
	return sqrt(n * (n + 1) - k * (k + 1));
}

double DiagSR(int N, int K, double Eaa, double Ebb, double Ecc, double J)
{
    //Hirota's with resimplification
    double val = C(J, N) / (2 * N * (N + 1)) * (Ecc * K * K + Ebb * (N * ( N + 1) - K * K));

    return val;
}

double NminusOneSR(int N, int K, double Eaa, double Ebb, double Ecc)
{
    //Hirota's with resimplification
    double val = K / (2 * N) * sqrt(N * N - K * K) * (Ebb - Ecc);

    return val;
}

double CentDist(double DN, double DNK, double DK, int N, int K)
{
    double val = -1.0 * DN * N * N * (N + 1) * (N + 1) - DNK * N * (N + 1) * K * K - DK * K * K * K * K;

    return val;
}

/* Model Name */
char* Name()
{
  return "NO3_Fourfold_beta";
}

/*-------------------------State types-------------------------*/
/* Number of state types. The actual number of the state types is 2:
the singlefold oblate symmetric top state (from Ming Wei Chen) and
the fourfold degenerate state*/
UINT StaNum()
{
  return 2;
}

char* StName( UINT nStateType )
{
	switch (nStateType)
	{
	case 0:
		return "Isolated State";
	case 1:
	default:
		return "Fourfold";
	}
}

/*--------------------------Constants--------------------------*/
/* Number of constants in a given state type. Constants are defined
above.*/
UINT ConNum( UINT nStateType )
{
	switch(nStateType)
	{
	case 0:
		return N_CONST_SINGLE;
	case 1:
		return N_CONST_TWOFOLD;
	default : return 0;
	}

}

/* Constant names*/
char* ConName( UINT nStateType, UINT nConst )
{
	switch(nStateType)
	{
		//FROM MWC
	case 0:
		switch(nConst)
		{
			case 0 : return "C";  //A
			case 1 : return "B";  //B
			case 2 : return "null";  //C
			case 3 : return "Dz";
			case 4 : return "Dx";
			case 5 : return "D3";
			case 6 : return "Dk";
			case 7 : return "Djk";
			case 8 : return "Dj";
			case 9 : return "SDk";
			case 10 : return "SDj";
			case 11 : return "Ecc";  //Eaa
			case 12 : return "Ebb";  //Ebb
			case 13 : return "null";  //Ecc
			case 14 : return "1/2(Eab+Eba)";
			case 15 : return "1/2(Ebc+Ecb)";
			case 16 : return "1/2(Eac+Eca)";
			case 17: return "planarity (0=off,1=C, 2=Dk, 3=both)";
			default : /* Non of the defined*/
			return "Constant";
		}//end switch nConst
	case 1:
		switch(nConst)
		{
		case 0 : return "Bzz:A1";
		case 1 : return "Byy:A1";
		case 2 : return "Bxx:A1";
		case 3 : return "Bzz:A2";
		case 4 : return "Byy:A2";
		case 5 : return "Bxx:A2";
		case 6 : return "Bzz:E";
		case 7 : return "Byy:E";
		case 8 : return "Bxx:E";
		case 9 : return "epsilon_1";
		case 10 : return "h1:A1_E";
		case 11 : return "h1:A2_E";
		case 12 : return "h1:E";
		case 13 : return "zeta_t:E";
		case 14 : return "zeta_t:A1A2";
		case 15 : return "a0_zeta_e_dz:E";
		case 16 : return "a_zeta_e_dz:A1A2";
		case 17 : return "a0_zeta_e_d";
		case 18 : return "BzzA1A2_zeta_t_A1A2";
		case 19 : return "Delta_E_A1_A2";
		case 20 : return "Delta_E_A1_E";
		case 21 : return "Epsilon_xxyy";//epsilon aa + bb
		case 22 : return "Epsilon_zz";//epsilon cc
		case 23 : return "DN";
		case 24 : return "DNK";
		case 25 : return "DK";
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
		pdConst[0] = 0.2286274;/* Ground State A*/
		pdConst[1] = 0.4585445;/* Ground State B*/
		pdConst[2] = 0.0;/* Ground State C*/
		pdConst[6] = 1.047e-6;
		pdConst[7] = -2.062e-6;
		pdConst[8] = 1.0880e-6;
		pdConst[11] = 0.00074;
		pdConst[12] = -0.01642;
		return;

	case 1 :
		pdConst[0] = 0.21;//Bzz  A1
		pdConst[1] = 0.43;//Byy  A1
		pdConst[2] = 0.43;//Bxx  A1
		pdConst[3] = 0.0;//Bzz  A2
		pdConst[4] = 0.0;//Byy  A2
		pdConst[5] = 0.0;//Bxx  A2
		pdConst[6] = 0.0;//Bzz  E
		pdConst[7] = 0.0;//Byy  E
		pdConst[8] = 0.0;//Bxx  E
		pdConst[21] = 0.015;
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
  return 1; /* J*/
}

/* Number of "BAD" Quantum Numbers */
UINT QNnumB(UINT nStateType)
{
	switch(nStateType)
	{
		case 0 : return 3;//N, K, P
		case 1 : return 4;//N, K, P, Symm
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
		case 3 : return "P";
		default : return "QN";/* Non of the defined*/
		}

	case 1 :
		switch(nQNumb)
		{
		case 0 : return "J";
		case 1 : return "N";
		case 2 : return "K";
		case 3 : return "P";
		case 4 : return "Symm";
		default : return "QN";
		}
	}
}

/* All Quantum Numbers are stored as integer, but we should know when
   they are actually half integer. In that case we store doubled value and
   base is "2" */
UINT QNBase( UINT nStateType, UINT nQN )
{
	switch(nStateType)
	{
	case 0 :
		switch (nQN)//j is half integer
		{
		case 0 : return MBASE;//For J
		case 3 : return MBASE;//for parity in MWC's model
		default : return 1;
		}
	case 1 :
		switch(nQN)
		{
		case 0 : return MBASE;//For J
		default : return 1;
		}
	}
}
/* Quantum number set that is first in the loop */
void StartQN( UINT nStateType, int* pnQNd )
{
  pnQNd[0]=1; /* This corresponds to J=1/2. */
}

/* The quantum number loop iteration (continue) condition. The Jmax value is specified in
as an integer, which is no big deal since it only sets the limit, but for comparison purpose
we need to convert it to doubled form.*/
BOOL ContQN( UINT nStateType, int* pnQNd, double* pdParam )
{
  return(pnQNd[0] <= 2*pdParam[7]);//changed from < to <= to match MWC, I think he's right
}

/* The quantum number loop iteration(next set in the loop). */
void NextQN( UINT nStateType, int* pnQNd, double* pdParam )
{
  pnQNd[0]+=2;
  return;
}

/* Because the fourfold model is fairly involved in terms of the quantum number
indexing. Therefore, to facilitate calculation of the matrix elements and inserting them into the matrix properly, we will
write two book-keeping functions. One translates the linear array index i or j into the
set of quantum numbers that the row or column corresponds to, and the other function does the opposite
thing, i.e. identifies i or j by the set of quantum numbers. By writing this into a single
function we avoid possible errors in the Hamiltonian function due to improper calculation of
the off-diagonal elements.*/

idx lIndex(int *pnQNd, int nx)
{
	int lN, hN, NK, spot;
	idx out;

	//two values of N for each J block, lN = lower and hN = higher value of N
	lN = (pnQNd[0]-1)/2;
	hN = (pnQNd[0]+1)/2;
	NK = 2 * (hN + lN + 1);//number of N and K combinations for a given symmetry sub-block
	int sBlock = nx / NK;
	if(sBlock < 2)//means it's an A1 or A2 level
	{
		out.P = 0;//Parity set to 0 for A1 and A2 levels
	}
	else
	{
		out.P = sBlock - 1;
	}
	out.Symm = sBlock;//this is to assign the symmetry labels. 0 for A1, 1 for A2, 2 for E-, 3 for E+
	spot = nx % NK;
	if(spot < 2 * hN + 1)
	{
		out.N = hN;//assigns N
		out.K = -1 * hN + spot;//assigns K
	}
	else
	{
		out.N = lN;//assigns N
		out.K = -1 * lN + spot - (2 * hN + 1);//assigns K
	}
	/* out is instance of the special structure with quantum numbers v, N, K.*/
	return out;
}


/*This function serves the same purpose as the lIndex function
and finds the quantum numbers for a given position in MWC's
isolated state Hamiltonian*/
idx MWCIndex(int *pnQNd, int nx)
{
	int lN, hN;
	idx out;
	//two values of N for each J block, lN = lower and hN = higher value of N
	lN = (pnQNd[0]-1)/2;
	hN = (pnQNd[0]+1)/2;
	out.Symm = 0;//since this is for isolated state the Symm doesn't matter
	out.P = 1;//same as Symm
	if(nx < 2 * hN + 1)
	{
		out.N = hN;
		out.K = -1 * hN + nx;
	}
	else
	{
		out.N = lN;
		out.K = -1 * lN + nx - (2 * hN + 1);
	}
	return out;
}

/*The second indexing function is for the reverse indexing.
Here, we don't check for the validity of quantum numbers, i.e.
nonnegative N, etc.*/
int lReverseIndex(int *pnQNd, int N, int K, int Symm)
{
  int NK = 2 * pnQNd[0] + 2;//NK is total number of N/K combinations possible for a given J value
  int nMark;//will be the index for where in a given symmetry block we are
  if(N > pnQNd[0]/2)//means N is the larger of two possible values for a given J
  {
	  nMark = K + N;
  }
  if(N == pnQNd[0]/2)//means N is the smaller of two possible values
  {
	  nMark = 3* N + 3 + K;//2(N + 1) + 1 + K + N
  }
  int nx = NK * Symm + nMark;
  return nx;
}

/*------------------------Main Part------------------------*/
/* Assign() function is not required for the model to operate but it helps to give meaningful
assignments. IMPORTANT NOTE: in this model assigment is given in terms of Wang components,
therefore K assumes only non-negative values. If such an assignment is implemented, DO NOT
use the MaxDltK option at all, since this will confuse the core and will result in
missing transitions. Also, in this program the symmetry labels will be given to the
eigenvectors, which currently are not used anywhere except in the user-end output.*/

/* NOTE 2: one of the issues is that the mechanism implemented in Specview deals with one
eigenvector at a time. In weakly coupled systems this won't be an issue, however it is
not clear how the assignment (labeling) function will behave in near-degenerate cases. */

void Assign( UINT nStateType, double* pdStWF, int* StQN, int StNum )
{
	if(nStateType == 0)//MWC's Assignment function
	{
		int N = (StQN[0]+1)/2;
		int indx = 0;
		double max = fabs(pdStWF[0]);
		for (int i = 1; i<4*N; i++)
		{
			if (fabs(pdStWF[i])>max)
			{
				max = fabs(pdStWF[i]);
				indx = i;
			}
		}

		if (indx<=2*N)
		{
			StQN[1] = N;
			StQN[2] = abs(indx-N);
			if (indx-N == 0)
			{
				if (PAR(N) > 0) StQN[3] = 0;
				else StQN[3] = 3;
				return;
			}
	//		if (fabs(max-fabs(pdStWF[2*N-indx]))<0.5*max)
	//		{
	//			if (pdStWF[2*N-indx]*pdStWF[indx]*PAR(indx)<0) StQN[3] = -1;
	//			else StQN[3] = 1;
				if (pdStWF[2*N-indx]*pdStWF[indx]>0)
				{
					if (PAR(StQN[2]) > 0)
					{
						if (PAR(N)>0) StQN[3] = 0;
						else StQN[3] = 3;
					}
					else
					{
						if (PAR(N)>0) StQN[3] = 1;
						else StQN[3] = 2;
					}
				}
				else
				{
					if (PAR(StQN[2]) > 0)
					{
						if (PAR(N)>0) StQN[3] = 3;
						else StQN[3] = 0;
					}
					else
					{
						if (PAR(N)>0) StQN[3] = 2;
						else StQN[3] = 1;
					}
				}

	//		}
	//		else StQN[3] = 0;
		}//end if indx<= 2 * N
		else
		{
			StQN[1] = N-1;
			StQN[2] = abs(indx - 3*N);
			if (indx - 3*N == 0)
			{
				if (PAR(N-1)>0) StQN[3] = 0;
				else StQN[3] = 3;
				return;
			}
	//		if (fabs(max-fabs(pdStWF[6*N-indx]))<0.5*max)
	//		{
	//			if (pdStWF[indx]*pdStWF[6*N-indx]*PAR(indx-1)<0) StQN[3] = -1;
	//		if (pdStWF[indx]*pdStWF[6*N-indx]>0) StQN[3] = -1;
			if (pdStWF[indx]*pdStWF[6*N-indx]>0)
			{
				if (PAR(StQN[2]) > 0)
				{
					if (PAR(N-1)>0) StQN[3] = 0;
					else StQN[3] = 3;
				}
				else
				{
					if (PAR(N-1)>0) StQN[3] = 1;
					else StQN[3] = 2;
				}
			}
			else
			{
				if (PAR(StQN[2]) > 0)
				{
					if (PAR(N-1)>0) StQN[3] = 3;
					else StQN[3] = 0;
				}
				else
				{
					if (PAR(N-1)>0) StQN[3] = 2;
					else StQN[3] = 1;
				}
			}
	//		}
	//		else StQN[3] = 0;
		}//end else
		//StQN[4] = 1;//just assigns symmetry to A2 for ground state
		return;
	}//MWC's assign

	if(nStateType == 1)//Working assignment function for Fourfold states
	{
	    int j = 0;
		idx QN;
		int id, DIM;
		DIM = 8*(StQN[0] + 1);
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
		StQN[1] = QN.N;
		StQN[2] = QN.K;
		StQN[3] = QN.P;
		StQN[4] = QN.Symm;
		return;
	}//My simple assign function
}//end assign

/* Get Hamiltonian matrix size for a given J block */
UINT HamSize( UINT nStateType, int* pnQNd )
{
  return (nStateType==1)? 8*(pnQNd[0] + 1): 2*(pnQNd[0] + 1);  /* number-of-ev-states*2*(2J+1), changed to 8* for nStateType = 1*/
}

/* Building Hamiltonian blocks for a given J, includes code for all supported states*/
void Hamilt( UINT nStateType,
	     double* pdConst,
	     int* pnQNd,
	     double** ppdH,int** ppnQNm)
{
	if(nStateType == 0)//MWC's Hamiltonian
	{
		double F( double J, double K );
  double g(double J, double K);

//Van Vleck notation:
  int i,j,Num, N, K;
  double J;//, P, U;
  double A,B,Cc,Dz,Dx,D3,Dk,Djk,Dj,SDk,SDj,Eaa,Ebb,Ecc,Eabba,Ebccb,Eacca,cons1,cons2,cons3;
  double planarity;
  B = pdConst[1];
  Cc = B;

  Dz = pdConst[3];
  Dx = pdConst[4];
  D3 = pdConst[5];
  Djk = pdConst[7];
  Dj = pdConst[8];
  SDk = pdConst[9];
  SDj = pdConst[10];
  Eaa = pdConst[11];
  Ebb = pdConst[12];
  Ecc = Ebb;
  Eabba = pdConst[14];
  Ebccb = pdConst[15];
  Eacca = pdConst[16];
  planarity = pdConst[17];

  if (planarity==1)
  {
	  A = B/2;/* constrain to planarity */
	  Dk = pdConst[6];
  }
  else if (planarity==2)
  {
	  A = pdConst[0];
	  Dk = -0.25 * ( 2 * Dj + 3 * Djk ); /* constrain to planarity */
  }
  else if (planarity==3)
  {
	  A = B/2;/* constrain to planarity */
	  Dk = -0.25 * ( 2 * Dj + 3 * Djk ); /* constrain to planarity */
  }
  else
  {
	  A = pdConst[0];
	  Dk = pdConst[6];
  }

  cons1 = (B+Cc)*0.5;//A = pdConst[0];
  cons2 = A-(B+Cc)*0.5;//B = pdConst[1];
  cons3 = (B-Cc)*0.5;//C = pdConst[2];

  double a0 = -(Eaa+Ecc+Ebb)/3;
  double a = -(2*Eaa-Ebb-Ecc)/6;
  double b = (Ecc-Ebb)/2;

  double c = -Ebccb;
  double d = -Eabba;
  double e = -Eacca;

	double Dcsk, Dcsnk, Dcskn, Dcsn, Dsn, Dsk;
	Dcsk = Dcsnk = Dcskn = Dcsn = Dsn = Dsk = 0;

	J = ((double)pnQNd[0])/2; /* J QUANTUM NUMBER*/
	Num = (pnQNd[0] + 1);

	for(i = 0; i < Num*2; i++)  /* INITIALIZE THE HAMILTONIAN MATRIX*/
		for(j = 0; j < Num*2;j++) ppdH[i][j] = 0.0;

	N = Num/2;
	for( i = 0; i < Num + 1; i++ )
	{
		K = i - N; /* K QUANTUM NUMBER*/

    /* DIAGONAL ELEMENTS*/
//restore
		ppdH[i][i] = cons2*K*K + cons1*N*(N+1) + Dcsk*K*K*K*K +
			 + Dz*K +
			D3*K*K*K - Dk*K*K*K*K - Djk*K*K*N*(N+1) - Dj*N*N*(N+1)*(N+1)
			-0.5*a0*C(J,N)+(3*K*K-N*(N+1))*a*Fe(J,N);//(1)
//restore
		if(i<Num)
		{//(2)
			ppdH[i][i+1] = -0.5*Dx*F(N, K)
				+(K+.5)*F(N,K)*Fe(J,N)*(d+e);
			ppdH[i+1][i] = ppdH[i][i+1];//(2)
		}
//restore
		if(i<Num-1)
		{//(3)
			ppdH[i+2][i] = ppdH[i][i+2] = cons3/2*F(N, K)*F(N, K+1) +
			- (0.5*SDk*(K*K+(K+2)*(K+2))+SDj*N*(N+1))*F(N, K)*F(N, K+1)
			+F(N, K)*F(N, K+1)/2*Fe(J,N)*(b+c);
		}
		if (abs(K)<N)
		{
//restore
			if (N>1) ppdH[i+Num][i+Num] = cons2*K*K + cons1*N*(N-1) +
				+ Dz*K +
			D3*K*K*K - Dk*K*K*K*K - Djk*K*K*N*(N-1) - Dj*N*N*(N-1)*(N-1)
			-0.5*a0*C(J,N-1)+(3*K*K-N*(N-1))*a*Fe(J,N-1);//(1)
			if (N == 1) ppdH[i+Num][i+Num] = 0;
//restore
			ppdH[i+Num][i] = ppdH[i][i+Num] =
				+3./2*K*sqrt(double (N*N-K*K))*a/N;//(4)
//-
			ppdH[i+Num][i+1] = ppdH[i+1][i+Num] =
				+(N-2*K-1)*g(N, -K-1)/4/N*(d-e);//(5)

			if(i<Num-1)
			{
				ppdH[i+Num][i+Num+1] = 0.5*Dx*F(N-1, K)
					+(K+.5)*F(N-1,K)*Fe(J,N-1)*(d+e);
				ppdH[i+Num+1][i+Num] = ppdH[i+Num][i+Num+1];//(2)

				ppdH[i][i+Num+1] = ppdH[i+Num+1][i] =
					+(N+2*K+1)*g(N,K)/4/N*(d+e);//(5)
//-
				ppdH[i+2][i+Num] = ppdH[i+Num][i+2] =
					-F(N, -K-2)/4*g(N, -K-1)/N*(b-c);//(6)
//-
			}

			if (abs(K+2)<N)
			{
//restore				//(3)
				ppdH[i+Num][i+Num+2] = ppdH[i+Num+2][i+Num] = cons3/2*F(N-1, K)*F(N-1, K+1)
				  - (0.5*SDk*(K*K+(K+2)*(K+2))+SDj*N*(N-1))*F(N-1, K)*F(N-1, K+1)
				  +F(N-1, K)*F(N-1, K+1)/2*Fe(J,N-1)*(b+c);
				//(6)
				ppdH[i+Num+2][i] = ppdH[i][i+Num+2] =
					+F(N,K)/4*g(N,K+1)/N*(b+c);//(6)
			}
		}
		ppnQNm[i][0] = N;
		ppnQNm[i][1] = K;
		ppnQNm[i][2] = 1;
		if (abs(K)<N)
		{
			ppnQNm[i+Num][0] = N-1;
			ppnQNm[i+Num][1] = K;
			ppnQNm[i+Num][2] = 1;
		}
	}
	}//end of nStateType == 0

	else if(nStateType == 1)//Fourfold Hamiltonian
	{
		double BzzE, BzzAO, BzzAT, ByyE, ByyAO, ByyAT, BxxE, BxxAO, BxxAT, epsilon, zetaEt, zetaAAt, hAOE,
			hATE, hE, azetaDE, azetaDAA, azetaD, bzzAOATzetat, dEAOAT, dEAOE, JJ, Exx, Eyy, Ezz, DN, DK, D_NK;

		int DIM, N, K, Par, Symm, dJ, NK, jj, coeff;
		idx QN;

		BzzAO = pdConst[0];
		ByyAO = pdConst[1];
		BxxAO = pdConst[2];
		BzzAT = pdConst[3];
		ByyAT = pdConst[4];
		BxxAT = pdConst[5];
		BzzE = pdConst[6];
		ByyE = pdConst[7];
		BxxE = pdConst[8];
		epsilon = pdConst[9];
		hAOE = pdConst[10];
		hATE = pdConst[11];
		hE = pdConst[12];
		zetaEt = pdConst[13];
		zetaAAt = pdConst[14];
		azetaDE = pdConst[15];
		azetaDAA = pdConst[16];
		azetaD = pdConst[17];
		bzzAOATzetat = pdConst[18];
		dEAOAT = pdConst[19];
		dEAOE = pdConst[20];
		Exx = pdConst[21];
		Eyy = pdConst[21];
		Ezz = pdConst[22];
		DN = pdConst[23];
		D_NK = pdConst[24];
		DK = pdConst[25];

		DIM = 8*(pnQNd[0] + 1);
		NK = DIM / 4;
		dJ = pnQNd[0];
		JJ = (double)pnQNd[0]/2.0;

		for(int ii = 0; ii < DIM; ii++)  /* INITIALIZE THE HAMILTONIAN MATRIX*/
		{
			for(int jjj = 0; jjj < DIM; jjj++)
			{
				ppdH[ii][jjj] = 0.0;
			}
		}

		for(int i = 0; i < DIM; i++)
		{
			//First, get Quantum Numbers
			QN = lIndex(pnQNd, i);

			K = QN.K;
			N = QN.N;
			Par = QN.P;
			Symm = QN.Symm;
			ppnQNm[i][0] = N;
			ppnQNm[i][1] = K;
			ppnQNm[i][2] = Par;
			ppnQNm[i][3] = Symm;

			//Vibronically Diagonal Elements First
			if(Symm == 0)//means its A1
			{
				ppdH[i][i] += BzzAO * K * K + 0.5 * (BxxAO + ByyAO) * (N * (N + 1) - K * K) + DiagSR(N, K, Eyy, Exx, Ezz, JJ) + CentDist(DN, D_NK, DK, N, K); //what I had originally for Hirota's formulas: DiagSR(N, K, Ezz, Exx, Eyy, JJ);//C.16
			}
			if(Symm == 1)//means its A2
			{
				ppdH[i][i] += BzzAT * K * K + 0.5 * (BxxAT + ByyAT) * (N * (N + 1) - K * K) + dEAOAT + DiagSR(N, K, Ezz, Exx, Eyy, JJ) + CentDist(DN, D_NK, DK, N, K);//C.16
			}
			if(Symm < 2 && dJ == 2 * N - 1)//means A1 or A2 and N is the larger of two values
			{
			    if(K < N && K > -1 * N)//means K isn't too large
                {
                    int o = lReverseIndex(pnQNd, N - 1, K, Symm);
                    ppdH[i][o] += NminusOneSR(N, K, Eyy, Exx, Ezz); //what I had originally for Hirota's formulas: NminusOneSR(N, K, Ezz, Exx, Eyy);
                    ppdH[o][i] = ppdH[i][o];
                }
			}
			if(Symm > 1)//means its E
			{
				if(dJ == 2 * N + 1)//means J = N + 1/2
				{
					ppdH[i][i] += BzzE * K * K + 0.5 * (BxxE + ByyE) * (N * (N + 1) - K * K) - (2 * BzzE * zetaEt -
						azetaDE / (dJ + 1)) * K + dEAOE + DiagSR(N, K, Eyy, Exx, Ezz, JJ) + CentDist(DN, D_NK, DK, N, K);//C.15a
				}
				else//for J = N - 1/2
				{
					ppdH[i][i] += BzzE * K * K + 0.5 * (BxxE + ByyE) * (N * (N + 1) - K * K) + (2 * BzzE * zetaEt -
						azetaDE / (dJ + 1)) * K + dEAOE + DiagSR(N, K, Eyy, Exx, Ezz, JJ) + CentDist(DN, D_NK, DK, N, K);//C.15a

					//can only have the matrix element with bra N - 1 here since it means ket's N is higher
					if(K < N && K > -1 * N)//first check to see if K value is not too large or small to be found in N - 1 set
					{
						jj = lReverseIndex(pnQNd, N - 1, K, Symm);
						ppdH[jj][i] += -1 * (azetaDE * sqrt((double)(N * N) - (double)(K * K))/(dJ + 1)) + NminusOneSR(N, K, Ezz, Exx, Eyy);//i + N - 1 should put me at same K for smaller N value ----C.15b
						ppdH[i][jj] = ppdH[jj][i];//since H is symmetric C.15b
					}
				}//end else

				if(-1 * K + 2 <= N)//make sure that K value for bra is not out of bounds
				{
					jj = lReverseIndex(pnQNd, N, -1 * K + 2, Symm);
					coeff = N - K + 2 + Symm - 1;//Symm - 1 is the parity value (1 for E- and 2 for E+)
					ppdH[jj][i] += hE * pow(-1.0, (double)coeff) * FNK(N, K - 1) * FNK(N, K - 2);//C.15c
					ppdH[i][jj] = ppdH[jj][i];//C.15c

					//now add the spin rotation terms to the same matrix elements
					int factor = 1;
					if(dJ == 2 * N - 1)
					{
						factor = -1;
					}

					ppdH[jj][i] += factor * coeff * epsilon / (dJ + 1) * FNK(N, K - 1) * FNK(N, K - 2);//C.18a
					ppdH[i][jj] = ppdH[jj][i];//C.18a
					if(dJ == 2 * N - 1 && -1 * K + 2 <= N - 1)//make sure that ket (column) is larger N
					{
						jj = lReverseIndex(pnQNd, N - 1, -1 * K + 2, Symm);
						ppdH[jj][i] += pow(-1.0, (double)coeff + 1.0) * epsilon / (dJ + 1) * sqrt((double)((N + K) * (N + K - 1) * (N + K - 2) * (N - K + 1)));//C.18b
						ppdH[i][jj] = ppdH[jj][i];//C.18b
					}
				}//end if
			}//end if Symm > 1 (meaning it's an E symm state)


			//Now onto the vibronically off diagonal terms
			//First A1 A2 terms
			if(Symm == 0)//means A1, take care of all terms here
			{
				jj = lReverseIndex(pnQNd, N, K, 1);
				ppdH[i][jj] += -2 * bzzAOATzetat * K;//C.19a
				ppdH[jj][i] = ppdH[i][jj];//C.19a

				if(dJ == 2 * N + 1)//means J = N + 1/2
				{
					ppdH[i][jj] += -1 * azetaDAA * K / (dJ + 1);//C.20a
					ppdH[jj][i] += -1 * azetaDAA * K / (dJ + 1);//C.20a

					//since J = N + 1/2 this is the smaller of the two N values, no need to do bounds check on K
					int jjjj = lReverseIndex(pnQNd, N + 1, K, 1);
					ppdH[jjjj][i] += azetaDAA * sqrt((JJ + 0.5) * (JJ + 0.5) - (double)(K * K)) / (double)(dJ + 1);//C.20b
					ppdH[i][jjjj] = ppdH[jjjj][i];//C.20b
				}
				else
				{
					ppdH[i][jj] += azetaDAA * K / (dJ + 1);//C.20a
					ppdH[jj][i] += azetaDAA * K / (dJ + 1);//C.20a

					if(K < N && K > -1 * N)
					{
						jj = lReverseIndex(pnQNd, N - 1, K, 1);
						ppdH[jj][i] += azetaDAA * sqrt((JJ + 0.5) * (JJ + 0.5) - (double)(K * K)) / (double)(dJ + 1);//C.20b
						ppdH[i][jj] = ppdH[jj][i];//C.20b
					}//end if
				}//end else
			}//end if Symm == 0

			//Next A1/A2 E terms
			if(Symm < 2)//means A1 or A2
			{
				if(K - 2 >= -1 * N && Symm == 0)//For A1 E
				{
					jj = lReverseIndex(pnQNd, N, K - 2, 2);//for E-
					ppdH[jj][i] += 1 / sqrt((double)2) * hAOE * FNK(N, K - 1) * FNK(N, K - 2);//C.22c
					ppdH[i][jj] = ppdH[jj][i];//C.22a

					jj = lReverseIndex(pnQNd, N, K - 2, 3);//for E+
					ppdH[jj][i] += 1 / sqrt((double)2) * hAOE * FNK(N, K - 1) * FNK(N, K - 2);//C.22c
					ppdH[i][jj] = ppdH[jj][i];//C.22a
				}
				if(K - 2 >= -1 * N && Symm == 1)//For A2 E
				{
					jj = lReverseIndex(pnQNd, N, K - 2, 2);//for E-
					ppdH[jj][i] += 1 / sqrt((double)2) * hATE * FNK(N, K - 1) * FNK(N, K - 2);//C.22c
					ppdH[i][jj] = ppdH[jj][i];//C.22a

					jj = lReverseIndex(pnQNd, N, K - 2, 3);//for E+
					ppdH[jj][i] += 1 / sqrt((double)2) * hATE * FNK(N, K - 1) * FNK(N, K - 2);//C.22c
					ppdH[i][jj] = ppdH[jj][i];//C.22a
				}

				if(-1 * K - 2 >= -1 * N && Symm == 0)
				{
					coeff = N + K + 2;
					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 2);//for A1 E-
					ppdH[jj][i] += pow(-1.0, (double)coeff + 1.0) / sqrt((double)2) * hAOE * FNK(N, K) * FNK(N, K + 1);//C.22d
					ppdH[i][jj] = ppdH[jj][i];//C.22d

					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 3);//for A1 E+
					ppdH[jj][i] += pow(-1.0, (double)coeff + 2.0) / sqrt((double)2) * hAOE * FNK(N, K) * FNK(N, K + 1);//C.22d
					ppdH[i][jj] = ppdH[jj][i];//C.22d
				}
				if(-1 * K - 2 >= -1 * N && Symm == 1)
				{
					coeff = N + K + 2;
					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 2);//for A2 E-
					ppdH[jj][i] += pow(-1.0, (double)coeff + 1.0 + (double)Symm) / sqrt((double)2) * hATE * FNK(N, K) * FNK(N, K + 1);//C.22d
					ppdH[i][jj] = ppdH[jj][i];//C.22d

					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 3);//for A2 E+
					ppdH[jj][i] += pow(-1.0, (double)coeff + 2.0 + (double)Symm) / sqrt((double)2) * hATE * FNK(N, K) * FNK(N, K + 1);//C.22d
					ppdH[i][jj] = ppdH[jj][i];//C.22d
				}
			}//end if Symm < 2

			if(Symm > 1)//means it's and E level
			{
				if(K + 2 <= N)
				{
					jj = lReverseIndex(pnQNd, N, K + 2, 0);//for A1 E
					ppdH[jj][i] += 1 / sqrt((double)2) * hAOE * FNK(N, K) * FNK(N, K + 1);//C.22a
					ppdH[i][jj] = ppdH[jj][i];//C.22a

					jj = lReverseIndex(pnQNd, N, K + 2, 1);//for A2 E
					ppdH[jj][i] += 1 / sqrt((double)2) * hATE * FNK(N, K) * FNK(N, K + 1);//C.22a
					ppdH[i][jj] = ppdH[jj][i];//C.22a
				}
				if(-1 * K - 2 >= -1 * N)
				{
					coeff = N - K + 2 + Par;
					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 0);//for A1 E
					ppdH[jj][i] += pow(-1.0, (double)coeff) / sqrt((double)2) * hAOE * FNK(N, K) * FNK(N, K + 1);//C.22b
					ppdH[i][jj] = ppdH[jj][i];//C.22b

					jj = lReverseIndex(pnQNd, N, -1 * K - 2, 1);//for A2 E
					ppdH[jj][i] += pow(-1.0, (double)coeff + 1.0) / sqrt((double)2) * hATE * FNK(N, K) * FNK(N, K + 1);//C.22b
					ppdH[i][jj] = ppdH[jj][i];//C.22b
				}
			}//end if Symm > 1
		}//end Hamiltonian for loop
	}//end if StateType == 1 (meaning its the E'' state)

	return;
}//end Hamilt



/* Get derivative of the Hamiltonian matrix on constant number k for a
   given good q.num. set */
void DHamilt( UINT nStateType,
              double* pdConst,
              int* pnQNd,
              double** ppdDH, int k )
{

	  UINT i, j;
	  double *pdIndicators, **ppdHbase, *pdOnes;
	  int **ppnQNm;
	  UINT nDim;
	  int lCN;
	  lCN = ConNum(nStateType);

	  pdIndicators = (double *)calloc(lCN , sizeof( double ) );
	  pdOnes = (double *)calloc(lCN , sizeof(double));

	  nDim = HamSize(nStateType,pnQNd);  /* Hamiltonian is (4*J+2)x(4*J+2) */

	  for ( int i = 0; i < lCN; i++ )
      {
          pdIndicators[i]=pdConst[i];
          pdOnes[i]=pdConst[i];
          //pdIndicators[i]=1.0;
          //pdOnes[i] = 1.0;
      }
	  //pdIndicators[k] +=1.0;
	  pdOnes[k] = 1.0;
	  pdIndicators[k] = 2.0;
	  ppnQNm = (int**)calloc( nDim, sizeof( int* ) );
	  ppdHbase = (double**)calloc(nDim, sizeof(double*));

	  for ( i = 0; i < nDim; i++ )
	  {
		ppnQNm[i] = (int*)calloc( QNnumB(nStateType), sizeof(int) );
		ppdHbase[i] = (double*)calloc( nDim, sizeof(double));
	  }
	//Hamilt (nStateType, pdConst, pnQNd, ppdHbase, ppnQNm);
	Hamilt (nStateType, pdOnes, pnQNd, ppdHbase, ppnQNm);
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
	free(pdOnes);
}

//equation C.39 in Dmitry's write up
double perpIntensity(int i, int* UpStQN, int* LoStQN, idx LoStIndex, int symm, int UpStateK, int UpStateN, double* pdLoStWF, double* pdUpStWF, double dWeightB)
{
      int j = lReverseIndex(UpStQN, UpStateN, UpStateK, symm);
      double coeff = pow(-1.0, (double)(UpStQN[0] + 1) / 2.0 + (double)UpStateK) / sqrt(2.0);
      double SixJ = TDM6J(UpStQN[0], UpStateN, LoStQN[0], LoStIndex.N);
      double ThreeJ_One = TDM3J(UpStateN, -1 * UpStateK, LoStIndex.N, LoStIndex.K, 1);
      double ThreeJ_Two = TDM3J(UpStateN, -1 * UpStateK, LoStIndex.N, -1 * LoStIndex.K, 1);
      double jjnn = JJNN(UpStQN[0], LoStQN[0], UpStateN, LoStIndex.N);
      double coeff_Two = pow(-1.0, LoStIndex.N - LoStIndex.K + 2.0 + 2.0 + symm - 1.0);//coefficient in front of second 3J symbol, p = symm - 1, k = j = 2.0
      double co = pow(-1.0, LoStIndex.N + UpStateN + 1);//takes care of 3J symbol being different from function
      double intensity = jjnn * coeff * SixJ * co * (ThreeJ_One + coeff_Two * ThreeJ_Two) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
      return intensity;
  }//end perpIntensity function

//eqtn C.32 in Dmitry's write up
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
	if(nLoStateType == 0 && nUpStateType == 0)//means both can use MWC's model
	{
		int f, nJL, nJU, nDelKmax/*, nKaL, nKaU,  SzL, SzU*/, NL, NU;
	double dNorm;
	double JL, JU;
	double static dWeightP, dWeightQ, dWeightR, dNormPQR, dWeightA, dWeightB, dWeightCsq, EMIN, RK;
	BOOL static bFlagA, bFlagBC;
	double NSSW(int N, int K);

  /* *pbCont = FALSE only for a first call*/
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





        dNormPQR  =  sqrt(pdParam[20]*pdParam[20]+pdParam[21]*pdParam[21]
			+pdParam[22]*pdParam[22]);

		dWeightP = fabs(pdParam[20])/dNormPQR;
		dWeightQ = fabs(pdParam[21])/dNormPQR;
		dWeightR = fabs(pdParam[22])/dNormPQR;

		dNorm = fabs(pdParam[9]) + fabs(pdParam[10]) + fabs(pdParam[11]);
		dNorm = (dNorm==0.0)?1.0:dNorm;
		dWeightA = sqrt(fabs(pdParam[9])/dNorm);     /* weight of a-type transition intensities.*/
		dWeightB = sqrt(fabs(pdParam[10]/dNorm));   /* weight of b-typetransition intensities.*/
		dWeightCsq = fabs(pdParam[11]/dNorm);    /* weight of c-type transition intensities.*/
		bFlagA =  (dWeightA != 0.0);        /*ladder.*/
		bFlagBC = ( dWeightB!=0.0 || dWeightCsq!=0.0);
	}
	*fInten = 0.0f;
	nJL = LoStQN[0];
	nJU = UpStQN[0];
//  nKaL = LoStQN[1];
//  nKaU = UpStQN[1];
//  SzU = UpStQN[3];
//  SzL = LoStQN[3];
	NL = (1+nJL)/2;
	NU = (1+nJU)/2;

	JL = 0.5*nJL;
	JU = 0.5*nJU;

	f = nJU - nJL;
	nDelKmax = (int)pdParam[8];

	if( abs(f)< 3 )
	{
    /* dJ <= 1, dK < dKmax  and Parities is diff.*/
		double  dTmp, SUM1, SUM2, SUM3, SUMAB, BF, LS, RINTE;
		int kL, KL;

		(*pnCount)++; /*  One more transition from     */
    /*  a current ground state level */

    /*  SUMMATIONS FOR a-type, b-type and c-type TRANSITIONS*/
		SUM1 = 0.0, SUM2 = 0.0, SUM3 = 0.0;
		double SJ00, SJ01, SJ10, SJ11, sqn, sqn1;
		sqn = sqrt(double (2*NL+1));
		sqn1 = -sqrt(double (2*NL-1));
		SJ00 = SJS(NL, JL, .5, JU, NU, 1)*sqn;
		SJ01 = SJS(NL, JL, .5, JU, NU-1, 1)*sqn;
		SJ10 = SJS(NL-1, JL, .5, JU, NU, 1)*sqn1;
		SJ11 = SJS(NL-1, JL, .5, JU, NU-1, 1)*sqn1;

		for ( KL=-NL; KL < NL+1; KL++)
		{
			kL = KL + NL;

			if (NU == NL)
			{
				if (bFlagA)
				{
					SUM1 += pdLoStWF[kL]*pdUpStWF[kL]*CG1(NL, 1., NU, KL, 0., KL)*
						SJ00
						*dWeightQ
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);
					if (abs(KL)<NL)
					{
						SUM1 += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU]*CG1(NL-1, 1., NU-1, KL, 0., KL)*
							SJ11
							*dWeightQ
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);
						SUM1 += pdLoStWF[kL+2*NL]*pdUpStWF[kL]*CG1(NL-1, 1., NU, KL, 0., KL)*
							SJ10
							*dWeightQ
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);;
						SUM1 += pdLoStWF[kL]*pdUpStWF[kL+2*NU]*CG1(NL, 1., NU-1, KL, 0., KL)*
							SJ01
							*dWeightQ
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);;
					}
				}

				if (bFlagBC && abs(KL-1)<=NL)
				{
					dTmp = pdLoStWF[kL]*pdUpStWF[kL-1]*CG1(NL, 1., NU, KL, -1., KL-1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);
					if (abs(KL)<NL)
					{
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL-1]*CG1(NL-1, 1., NU, KL, -1., KL-1.)*
							SJ10
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);

						if (abs(KL-1)<NL)
							dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU-1]*CG1(NL-1, 1., NU-1, KL, -1., KL-1.)*
								SJ11
								* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						        * ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);
					}

					if (abs(KL-1)<NL)
						dTmp += pdLoStWF[kL]*pdUpStWF[kL-1+2*NU]*CG1(NL, 1., NU-1, KL, -1., KL-1.)*
							SJ01
						    * ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);

					SUM2 += dTmp*dWeightQ;
					SUM3 += dTmp*dWeightQ;
				}


				if (bFlagBC && abs(KL+1)<=NL)
				{
					dTmp = pdLoStWF[kL]*pdUpStWF[kL+1]*CG1(NL, 1., NU, KL, 1., KL+1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);
					if (abs(KL)<NL)
					{
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+1]*CG1(NL-1, 1., NU, KL, 1., KL+1.)*
							SJ10
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);

						if (abs(KL+1)<NL)
							dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU+1]*CG1(NL-1, 1., NU-1, KL, 1., KL+1.)*
								SJ11
								* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						        * ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);
					}

					if (abs(KL+1)<NL)
						dTmp += pdLoStWF[kL]*pdUpStWF[kL+1+2*NU]*CG1(NL, 1., NU-1, KL, 1., KL+1.)*
							SJ01
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						    * ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);

					SUM2 -= dTmp*dWeightQ;
					SUM3 += dTmp*dWeightQ;
				}

			}

			if (NU+1 == NL)
			{
				if (bFlagA)
				{
					if (abs(KL)<NL)
					{
						SUM1 += pdLoStWF[kL]*pdUpStWF[kL-1]*CG1(NL, 1., NU, KL, 0., KL)*
							SJ00
							*dWeightP
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);

						SUM1 += pdLoStWF[kL+2*NL]*pdUpStWF[kL-1]*CG1(NL-1, 1., NU, KL, 0., KL)*
							SJ10
							*dWeightP
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);

						if (abs(KL)<NL-1)
							SUM1 += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU-1]*CG1(NL-1, 1., NU-1, KL, 0., KL)*
								SJ11
								*dWeightP
								* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);
					}
				}

				if (bFlagBC && abs(KL-1)<NL)
				{
					dTmp = pdLoStWF[kL]*pdUpStWF[kL-2]*CG1(NL, 1., NU, KL, -1., KL-1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);

					if (abs(KL)<NL)
					{
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL-2]*CG1(NL-1, 1., NU, KL, -1., KL-1.)*
							SJ10
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);

						if (abs(KL-1)<NL-1)
							dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU-2]*CG1(NL-1, 1., NU-1, KL, -1., KL-1.)*
								SJ11
								* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);
					}
					SUM2 += dTmp*dWeightP;
					SUM3 += dTmp*dWeightP;
				}



				if (bFlagBC && abs(KL+1)<NL)
				{
					dTmp = pdLoStWF[kL]*pdUpStWF[kL]*CG1(NL, 1., NU, KL, 1., KL+1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);
					if (abs(KL)<NL)
					{
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL]*CG1(NL-1, 1., NU, KL, 1., KL+1.)*
							SJ10
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);;

						if (abs(KL+1)<NL-1)
							dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU]*CG1(NL-1, 1., NU-1, KL, 1., KL+1.)*
								SJ11
								* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);;
					}
					SUM2 -= dTmp*dWeightP;
					SUM3 += dTmp*dWeightP;
				}

			}

			if (NU == NL+1)
			{
				if (bFlagA)
				{
					SUM1 += pdLoStWF[kL]*pdUpStWF[kL+1]*CG1(NL, 1., NU, KL, 0., KL)*
						SJ00*dWeightR
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);

					SUM1 += pdLoStWF[kL]*pdUpStWF[kL+1+2*NU]*CG1(NL, 1., NU-1, KL, 0., KL)*
						SJ01*dWeightR
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);

					if (abs(KL)<NL)
					{
						SUM1 += pdLoStWF[kL+2*NL]*pdUpStWF[kL+1+2*NU]*CG1(NL-1, 1., NU-1, KL, 0., KL)*
							SJ11*dWeightR
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL)==(int)pdParam[13]);
					}
				}

				if (bFlagBC)
				{
					dTmp = pdLoStWF[kL]*pdUpStWF[kL]*CG1(NL, 1., NU, KL, -1., KL-1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);
					if (abs(KL)<NL)
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU]*CG1(NL-1, 1., NU-1, KL, -1., KL-1.)*
							SJ11
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);
					if (abs(KL-1)<=NL)
						dTmp += pdLoStWF[kL]*pdUpStWF[kL+2*NU]*CG1(NL, 1., NU-1, KL, -1., KL-1.)*
							SJ01
							* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL-1)==(int)pdParam[13]);

					SUM2 += dTmp*dWeightR;
					SUM3 += dTmp*dWeightR;

					dTmp = pdLoStWF[kL]*pdUpStWF[kL+2]*CG1(NL, 1., NU, KL, 1., KL+1.)*
						SJ00
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);
					if (abs(KL)<NL)
						dTmp += pdLoStWF[kL+2*NL]*pdUpStWF[kL+2*NU+2]*CG1(NL-1, 1., NU-1, KL, 1., KL+1.)*
							SJ11
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);


					if (abs(KL+1)<=NL)
						dTmp += pdLoStWF[kL]*pdUpStWF[kL+2+2*NU]*CG1(NL, 1., NU-1, KL, 1., KL+1.)*
							SJ01
						* ((int)pdParam[14]==0? 1 : (KL)==(int)pdParam[12])
						* ((int)pdParam[15]==0? 1 : (KL+1)==(int)pdParam[13]);

					SUM2 -= dTmp*dWeightR;
					SUM3 += dTmp*dWeightR;

				}
			}
		}








//    SUMAB = dWeightA*SUM1 + 0.707*dWeightB*SUM3;
		SUMAB = dWeightA*SUM1 + 0.707*dWeightB*SUM2;


		LS = (2*JL+1)*(2*JU+1)*((SUMAB * SUMAB) + 0.5*dWeightCsq*SUM3*SUM3);

		BF = exp(-(dLoStEnerg - EMIN)/(RK*pdParam[0]))-exp(-(dUpStEnerg - EMIN)/(RK*pdParam[0]));   /* Boltzmann factor */
		RINTE = ((int)pdParam[23]==0? 1 : NSSW(LoStQN[1], LoStQN[2]))* BF * LS;                     /* Intensity */


		if(RINTE >= pdParam[3] && (LoStQN[1]==(int)pdParam[16] || !(int)pdParam[18]) &&
			(UpStQN[1]==(int)pdParam[17] || !(int)pdParam[19]))
		*fInten = (float)RINTE;
	}
    /* if( *pnCount >= 6) *pbCont = FALSE;*/

//  if (nJL==3 && nJU==3) *fInten = 1.0f;
//  else *fInten = 0.0f;
	return;
	}//end MWC's model intensity function

	if(nLoStateType == 0 && nUpStateType == 1)//means transition from isolated ground state to fourfold excited state
	{
		double dNorm, sum, BF, LS;
		int f, dJL, dJU, DIM, uNL, uNH, j, coeff;
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
		  DIM = 2*(LoStQN[0] + 1);//sets the bounds for the for loop over the GS EV components
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
				  LoStIndex = MWCIndex(LoStQN, i);
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

		  if(bFlagBC)//C.39
		  {
			  for(int symm = 2; symm < 4; symm++)
			  {
				  for(int i = 0; i < DIM; i++)
				  {
				      LoStIndex = MWCIndex(LoStQN, i);
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
			  }//end loop over E+/E-
		  }//end if bFlagBC

		  LS = sum*sum; /* "Line strength" */
		  BF = exp(-(dLoStEnerg - EMIN)/(RK*pdParam[0]))-exp(-(dUpStEnerg - EMIN)/(RK*pdParam[0]));   /* Boltzmann factor */
		  RINTE = ((int)pdParam[23]==0? 1 : NSSW(LoStQN[1], LoStQN[2])) * BF * LS;
		  *(pnCount)++;
		  if(RINTE >= pdParam[3])
		  {
			  *fInten = (float)RINTE;
		  }
		  return;
	}//end if loState == onefold and excited state == fourfold

    if(nLoStateType == 1 && nUpStateType == 1)//means both states are fourfold
    {
		double dNorm, sum, BF, LS;
		int f, dJL, dJU, DIM, uNL, uNH, j, coeff;
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
			  for(int i = DIM / 4; i < DIM / 2; i++)//only loop over A2 levels since that's the ground state symmetry
			  {
				  LoStIndex = lIndex(LoStQN, i);// MWCIndex(LoStQN, i);
				  //Now check for 3 possible components that the GS component can link with first checking to see if the bad QN's exist in the excited state and only check for A1 (Symm = 0)
				  coeff = LoStIndex.K + (UpStQN[0] + 1) / 2 + LoStIndex.N * 2;
				  if(LoStIndex.N - 1 == uNL || LoStIndex.N - 1 == uNH)
				  {
				      upStateN = LoStIndex.N - 1;
				      j = lReverseIndex(UpStQN, upStateN, LoStIndex.K, 0);
                      sum += pow(-1.0, coeff) * JJNN(UpStQN[0], LoStQN[0], upStateN, LoStIndex.N) * TDM6J(UpStQN[0], upStateN, LoStQN[0], LoStIndex.N)
                        * TDM3J(upStateN, LoStIndex.K, LoStIndex.N, LoStIndex.K, 0) * pdLoStWF[i] * pdUpStWF[j] * dWeightA;
                  }
				  if(LoStIndex.N == uNL || LoStIndex.N == uNH && LoStIndex.K != 0 && LoStIndex.N != 0)
				  {
				      upStateN = LoStIndex.N;
					  coeff += 1;
					  j = lReverseIndex(UpStQN, upStateN, LoStIndex.K, 0);
                      sum += pow(-1.0, coeff) * JJNN(UpStQN[0], LoStQN[0], upStateN, LoStIndex.N) * TDM6J(UpStQN[0], upStateN, LoStQN[0], LoStIndex.N)
                        * TDM3J(upStateN, LoStIndex.K, LoStIndex.N, LoStIndex.K, 0) * pdLoStWF[i] * pdUpStWF[j] * dWeightA;
                  }
				  if(LoStIndex.N + 1 == uNL || LoStIndex.N + 1 == uNH)
				  {
				      upStateN = LoStIndex.N + 1;
					  coeff += 2;
					  j = lReverseIndex(UpStQN, upStateN, LoStIndex.K, 0);
                      sum += pow(-1.0, coeff) * JJNN(UpStQN[0], LoStQN[0], upStateN, LoStIndex.N) * TDM6J(UpStQN[0], upStateN, LoStQN[0], LoStIndex.N)
                        * TDM3J(upStateN, LoStIndex.K, LoStIndex.N, LoStIndex.K, 0) * pdLoStWF[i] * pdUpStWF[j] * dWeightA;
                  }
			  }//end for loop
		  }//end if bFlagA == true

		  if(bFlagBC)//C.39
		  {
			  int co = 0;
			  for(int symm = 2; symm < 4; symm++)
			  {
				  for(int i = DIM / 4; i < DIM / 2; i++)//only loop over A2 levels
				  {
				      LoStIndex = lIndex(LoStQN, i);//MWCIndex(LoStQN, i);
				      co = pow(-1, LoStIndex.N * 2 + 1);//this is the coefficient for the 3J symbols
					  if(LoStIndex.N - 1 == uNL || LoStIndex.N - 1 == uNH)
					  {
					      co *= -1;
						  if(LoStIndex.K - 1 >= LoStIndex.N * -1)
						  {
						      upStateK = LoStIndex.K - 1;
							  j = lReverseIndex(UpStQN, LoStIndex.N - 1, LoStIndex.K - 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K - 1;//coefficient in front of the 6J symbol
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N - 1) * TDM6J(UpStQN[0], LoStIndex.N - 1, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N - 1, -1 * (LoStIndex.K - 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)//coefficient in front of second 3J symbol, p = symm - 1, k = j = 2.0
								  * TDM3J(LoStIndex.N - 1, -1 * (LoStIndex.K - 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
						  if(LoStIndex.K + 1 <= LoStIndex.N)
						  {
							  j = lReverseIndex(UpStQN, LoStIndex.N - 1, LoStIndex.K + 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K + 1;
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N - 1) * TDM6J(UpStQN[0], LoStIndex.N - 1, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N - 1, -1 * (LoStIndex.K + 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)
								  * TDM3J(LoStIndex.N - 1, -1 * (LoStIndex.K + 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
					  }

					  if(LoStIndex.N == uNL || LoStIndex.N == uNH && LoStIndex.N != 0)
					  {
						  if(LoStIndex.K - 1 >= LoStIndex.N * -1)
						  {
							  j = lReverseIndex(UpStQN, LoStIndex.N, LoStIndex.K - 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K - 1;
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N) * TDM6J(UpStQN[0], LoStIndex.N, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N, -1 * (LoStIndex.K - 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)
								  * TDM3J(LoStIndex.N, -1 * (LoStIndex.K - 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
						  if(LoStIndex.K + 1 <= LoStIndex.N)
						  {
							  j = lReverseIndex(UpStQN, LoStIndex.N, LoStIndex.K + 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K + 1;
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N) * TDM6J(UpStQN[0], LoStIndex.N, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N, -1 * (LoStIndex.K + 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)
								  * TDM3J(LoStIndex.N, -1 * (LoStIndex.K + 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
					  }

					  if(LoStIndex.N + 1 == uNL || LoStIndex.N + 1 == uNH)
					  {
					      co *= -1;
						  if(LoStIndex.K - 1 >= LoStIndex.N * -1)
						  {
							  j = lReverseIndex(UpStQN, LoStIndex.N + 1, LoStIndex.K - 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K - 1;
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N + 1) * TDM6J(UpStQN[0], LoStIndex.N + 1, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N + 1, -1 * (LoStIndex.K - 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)//coefficient in front of second 3J symbol, p = symm - 1, k = j = 2.0
								  * TDM3J(LoStIndex.N + 1, -1 * (LoStIndex.K - 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
						  if(LoStIndex.K + 1 <= LoStIndex.N)
						  {
							  j = lReverseIndex(UpStQN, LoStIndex.N + 1, LoStIndex.K + 1, symm);
							  coeff = (UpStQN[0] + 1) / 2 + LoStIndex.K + 1;
							  sum += pow(-1.0, coeff) / sqrt(2.0) * JJNN(UpStQN[0], LoStQN[0], LoStIndex.N, LoStIndex.N + 1) * TDM6J(UpStQN[0], LoStIndex.N + 1, LoStQN[0], LoStIndex.N)
								  * co * (TDM3J(LoStIndex.N + 1, -1 * (LoStIndex.K + 1), LoStIndex.N, LoStIndex.K, 1) + pow(-1.0, LoStIndex.N - LoStIndex.K + symm - 1 + 2.0 + 2.0)//coefficient in front of second 3J symbol, p = symm - 1, k = j = 2.0
								  * TDM3J(LoStIndex.N + 1, -1 * (LoStIndex.K + 1), LoStIndex.N, -1 * LoStIndex.K, 1)) * pdLoStWF[i] * pdUpStWF[j] * dWeightB;
						  }
					  }//end if N+1 = uNL/uNH

				  }//end loop over lostate wf
			  }//end loop over E+/E-
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
    }//end upstate == lowstate == 1

  }//end Inten function
  //}
