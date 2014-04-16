#ifndef __MODEL_H__
#define __MODEL_H__

#ifdef __cplusplus
extern "C" {
#endif

char* Name();   		 // Current Model Name

void InitCo( UINT nStateType, double *pdConst );
void InitPa( double *pdParam);

char*  QNName( UINT nStateType , UINT nQNumb ); // name of Quantum Number
char*  StName( UINT nStateType );               // name of State(Sig,Pi,..)
char*  ConName( UINT nStateType, UINT nConst ); // name of Constants
char*  ParName( UINT nParam );                  // name of Parameters

UINT  QNBase( UINT nStateType, UINT nQN ); // Number to devide QNumb.(1 or 2)

UINT  StaNum();                  // number of different state types
UINT  ConNum( UINT nStateType ); // number of Constants 
UINT  ParNum();                  // number of Parameters
     
UINT  QNnumG(UINT nStateType);   // number of "good" Quantum Numbers
UINT  QNnumB(UINT nStateType);   // number of "bad" Quantum Numbers

void  StartQN( UINT nStateType, int* pnQNd ); // initiall "good" QN set
BOOL  ContQN( UINT nStateType, int* pnQNd, double* pdParam );//Cont. or stop
void  NextQN( UINT nStateType, int* pnQNd, double* pdParam ); // next "good" QN set

UINT  HamSize( UINT nStateType, int* pnQNd ); // Dimen. of H 

void  Assign( UINT nStateType, double* pdStWF, int* StQN, int StNum );

void  Hamilt( UINT nStateType, // e.g. 1 <=> Sigm, 2 <=> Pi, 3 <=> Delt.
	   double* pdConst, // array of  molecular Constants
	   int* pnQNd,      // set of "diagonal" Quantum Numbers
	   double** ppdH,int** ppnQNm); // Hamiltonian and "mixed" QN.
			   
void  DHamilt( UINT nStateType, // e.g. 1 <=> Sigm, 2 <=> Pi, 3 <=> Delt.
	double* pdConst, // array of  molecular Constants
	int* pnQNd, 	 // set of "diagonal" Quantum Numbers
	double** ppdDH, int k );// Hamiltonian deriv.<=> dH[i,j]/dConst[k]
				
void  Inten( int* pnCount, // Counter of alowed transition.
	  BOOL* pbCont, // FALSE if there is no more tran. from lo. State
	  double* pdParam,   // array of  Parameters
	  double dLoStEnerg, // energy of the lower state
	  UINT nLoStateType,
	  double* pdLoStWF, int* LoStQN, // wavefunc. & assigned quantum
	  double dUpStEnerg,
	  UINT nUpStateType,
	  double* pdUpStWF, int* UpStQN, // numbers of lo. and up. states
	  float* fInten );   // Intencity of transition( 0 if prohibited)

#ifdef __cplusplus
}
#endif

#endif //__MODEL_H__
