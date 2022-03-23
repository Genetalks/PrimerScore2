#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <float.h> /* ! mul ei ole float.h-d includes DBL_MAX */

#if defined(__sun)
#include <ieeefp.h>
#endif

#include "alignment_TM.h"
#include "thal_parameters.h"

/* table where bp-s enthalpies, that retrieve to the most stable Tm, are saved */
#ifdef EnthalpyDPT
# undef EnthalpyDPT
#endif
#define EnthalpyDPT(i, j) enthalpyDPT[(j) + ((i-1)*len3) - (1)]

/* table where bp-s entropies, that retrieve to the most stable Tm, are saved */
#ifdef EntropyDPT
# undef EntropyDPT
#endif
#define EntropyDPT(i, j) entropyDPT[(j) + ((i-1)*len3) - (1)]

#define MAX_LOOP  30

static void reverse(unsigned char *s);
static double Ss(int i, int j, int k); /* returns stack entropy */
static double Hs(int i, int j, int k); /* returns stack enthalpy */
void alignment_tm(const unsigned char* align0, const unsigned char* seq1, const unsigned char* seq2, const thal_args* a, thal_results* o, double *tm, double *dG, double *dS, double *dH);
void SH_Ldange(int i, int j, double* EntropyEnthalpy);
static void SH_Rdange(int i, int j, double* EntropyEnthalpy);
static int length_unsig_char(const unsigned char * str); /* returns length of unsigned char; to avoid warnings while compiling */
static void RSH(int i, int j, double* EntropyEnthalpy);
char * usage_out();

static double* enthalpyDPT; /* matrix for values of enthalpy */
static double* entropyDPT; /* matrix for values of entropy */

const char *path;
char *endptr; /* reading input */
int i; /* index */
static int len1, len2, lena, len3;
static unsigned char *numSeq1, *numSeq2, *align; /* same as seq1 and seq2 but converted to numbers */

int maxLoop = MAX_LOOP;
void caculate_alignment_TM(const unsigned char* align0, const unsigned char* seq1, const unsigned char* seq2, double mv, double dv, double dntp, double dna_conc, double temp, double *tm, double *dG, double *dS, double *dH)
{   
   thal_args a;
   thal_results o;
   o.sec_struct=NULL;
   set_thal_default_args(&a);

    a.mv = mv;
    a.dv = dv;
    a.dna_conc = dna_conc;
    a.temp = temp + ABSOLUTE_ZERO;
    a.dntp = dntp;
    a.maxLoop = maxLoop;
            
	/* END reading INPUT */
   /* read default thermodynamic parameters */
   thal_parameters thermodynamic_parameters;
   thal_set_null_parameters(&thermodynamic_parameters);

   thal_load_parameters(path, &thermodynamic_parameters, &o);
   get_thermodynamic_values(&thermodynamic_parameters, &o);

   lena = length_unsig_char(align0);
   len1 = length_unsig_char(seq1);
   len2 = length_unsig_char(seq2);
   len3 = len2;
   if(len1 != lena || len2 !=lena){
	   fprintf(stderr, "Wrong: lengths of %s, %s or %s are not the same!\n", align0, seq1, seq2);
	   exit(-1);
   }
	alignment_tm(align0, seq1, seq2, &a, &o, tm, dG, dS, dH);
	
	/* cleanup */
   destroy_thal_structures();
   thal_free_parameters(&thermodynamic_parameters);

}


void alignment_tm(const unsigned char* align0, const unsigned char* seq1, const unsigned char* seq2, const thal_args* a, thal_results* o, double *tm, double *dG, double *dS, double *dH){
	int i, j;
	int len;
	double H, S, G, t;
    double* SH;
	int end1, end2;
    len = length_unsig_char(align0);
    dplx_init_H = 200;
    dplx_init_S = -5.7;
	RC = R  * log(a->dna_conc/4000000000.0);

   /* convert nucleotides to numbers */
   numSeq1 = (unsigned char*) safe_realloc(numSeq1, len + 2, o);
   numSeq2 = (unsigned char*) safe_realloc(numSeq2, len + 2, o);
   align = (unsigned char*) safe_realloc(align, len + 2, o);

	/*** Calc part of the salt correction ***/
   saltCorrection=saltCorrectS(a->mv,a->dv,a->dntp); /* salt correction for entropy, must be multiplied with N, which is
                                                   the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */

   reverse(seq2); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
   for(i = 1; i <= len; ++i) {
	   if(align0[i-1]=='#'){ 
		   numSeq1[i] = 0; // tm of 0+1(numSeq1=0, numSeq2=1) is 1 less than tm of 0+0;
		   numSeq2[i] = 1; // tm of 1+1(numSeq1=1, numSeq2=1) is 2-3 less than tm of 0+0;
	   }else{
		   /* A-0, C-1, G-2, T-3 */
		   numSeq1[i] = str2int(seq1[i - 1]); 
		   numSeq2[i] = str2int(seq2[i - 1]);
	   }
	   align[i] = align0[i-1];
   }
   numSeq1[0] = numSeq1[len + 1] = 4; /* mark as end */
   numSeq2[0] = numSeq2[len + 1] = 4; /* mark as N-s */
   align[0] = align[len + 1] = '\0'; 
	
	/*** caculate tm ***/
	S = dplx_init_S;
	H = dplx_init_H;
    SH = (double*) safe_malloc(2 * sizeof(double), o);
	int iterminal = -1;
	int N=0;
	for (i = 1; i <= len; i++) {
		if(align[i] != '|'){
			continue;
		}
		N++;
		if(iterminal==-1){ /* Terminal AT penalty */
			S += atPenaltyS(numSeq1[i], numSeq2[i]);
			H += atPenaltyH(numSeq1[i], numSeq2[i]);
		}
		iterminal = i;
	
		/* left */
		if(align[i-1] == '#'){ /* left dange */
			SH_Ldange(i, i, SH);
			S += SH[0];
			H += SH[1];
//		}else if(align[i-1] == '*'){ /* left mismatch */
//			S += tstack2Entropies[numSeq2[i]][numSeq2[i-1]][numSeq1[i]][numSeq1[i-1]];
//			H += tstack2Enthalpies[numSeq2[i]][numSeq2[i-1]][numSeq1[i]][numSeq1[i-1]];
		}else if(align[i-1] == '|'){ /* left align */
			S += stackEntropies[numSeq1[i-1]][numSeq1[i]][numSeq2[i-1]][numSeq2[i]];
			H += stackEnthalpies[numSeq1[i-1]][numSeq1[i]][numSeq2[i-1]][numSeq2[i]];
		}
	
		if(align[i+1] == '#'){ /* right dange */
			SH_Rdange(i, i, SH);
			S += SH[0];
			H += SH[1];
//		}else if(align[i+1]== '^' || align[i+1] =='-' || (align[i+1]=='*' && i+2<len && align[i+2]=='*')){ /* loop */
		}else if(align[i+1]== '^' || align[i+1] =='-' || align[i+1] == '*'){ /* loop */
			end1 = end2 = i;
			int off1=0;
			int off2=0;
			for(j=i+1; j<len+1; j++){
				if(align[j] =='^'){
					end1++;
					off2++;
				}else if(align[j] == '-'){
					off1++;
					end2++;
				}else if(align[j] == '*'){
					end1++;
					end2++;
				}else if(align[j] == '|'){
					end1++;
					end2++;
					break;
				}
			}
			calc_bulge_internal(i, i, end1, end2, SH, 0, a->maxLoop, off1, off2);
			i=j-1;
			S += SH[0];
			H += SH[1];
//		}else if(align[i+1]=='*'){ /* right mismatch */
//			S += stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[i]][numSeq2[i+1]];
//			H += stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[i]][numSeq2[i+1]];
		}
	//	printf("i=%d, dS=%f, dH=%f, %c\n", i, S, H, align[i]); 
	}
	S += atPenaltyS(numSeq1[iterminal], numSeq2[iterminal]);
	H += atPenaltyH(numSeq1[iterminal], numSeq2[iterminal]);

	/* caculate G and T */
	N--;
	t = ((H) / (S + (N * saltCorrection) + RC)) - ABSOLUTE_ZERO;
    G = (H) - ((ABSOLUTE_ZERO+37) * (S + (N * saltCorrection)));
    S = S + (N * saltCorrection);
	
	//printf("dS=%f, dH=%f, dG=%f, t=%f\n", S, H, G, t);

	*dS = S;
	*dH = H;
	*dG = G;
	*tm = t;
}


void SH_Ldange(int i, int j, double* EntropyEnthalpy){
   double S2, H2, T2, G2;
   S2 = -1.0;
   H2 = -_INFINITY;
   T2 = -_INFINITY;

   /** If there is two dangling ends at the same end of duplex **/
   if((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1 ) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]]) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]] +
        dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
   //   G2 = H2 - TEMP_KELVIN*S2;
   //   if(!isFinite(H2) || G2>0) {
   //      H2 = _INFINITY;
   //      S2 = -1.0;
   //  G2 = 1.0;
   //   }
   //   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq2[j]][numSeq2[j - 1]][numSeq1[i]];
      G2 = H2 - TEMP_KELVIN*S2;
   //   if(!isFinite(H2) || G2>0) {
   //      H2 = _INFINITY;
   //      S2 = -1.0;
   //      G2 = 1.0;
   //   }
    //  T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   } else if ((bpIndx(numSeq1[i-1], numSeq2[j-1]) != 1) && isFinite(dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq2[j]][numSeq1[i]][numSeq1[i - 1]];
   //    G2 = H2 - TEMP_KELVIN*S2;
   //   if(!isFinite(H2) || G2>0) {
   //      H2 = _INFINITY;
   //      S2 = -1.0;
   //      G2 = 1.0;     
   //   }
   //   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   }
   EntropyEnthalpy[0] = S2;
   EntropyEnthalpy[1] = H2;
	return;
}

static void SH_Rdange(int i, int j, double* EntropyEnthalpy){
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   S1 = H1 = H2 = _INFINITY;
   H1 = T1 = T2 = -_INFINITY;
/* two dangle: */
  /*     / */ 
  /* ---/ */ 
  /* ---\ */ 
  /*     \ */ 
   if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    //G2 = H2 - TEMP_KELVIN*S2;
    //  if(!isFinite(H2) || G2>0) {
    //     H2 = _INFINITY;
    //     S2 = -1.0;
    // G2 = 1.0;
    //  }
    //  T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   }
/* one dangle: */
  /*     / */ 
  /* ---/ */ 
  /* --- */ 
  /*      */ 
   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
    //  G2 = H2 - TEMP_KELVIN*S2;
    //  if(!isFinite(H2) || G2 >0) {
    //     H2 = _INFINITY;
    //     S2 = -1.0;
	//	 G2 = 1.0;
    //  }
    //  T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   }
/* one dangle: */
  /*      */ 
  /* --- */ 
  /* ---\ */ 
  /*     \ */ 
   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    //  G2 = H2 - TEMP_KELVIN*S2;
    //  if(!isFinite(H2) || G2>0) {
    //     H2 = _INFINITY;
    //     S2 = -1.0;
	//	 G2 = 1.0;
    //  }
    //  T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
    }

    EntropyEnthalpy[0] = S2;
	EntropyEnthalpy[1] = H2;
    return;
}

static void 
calc_bulge_internal(int i, int j, int ii, int jj, double* EntropyEnthalpy, int traceback, int maxLoop, int off1, int off2)
{
   int loopSize1, loopSize2, loopSize;
   double S,H,G1,G2;
   int N, N_loop;
   double* SH;
   SH = (double*) safe_malloc(2 * sizeof(double), 0);
   SH[0] = -1.0;
   SH[1] = _INFINITY;
   S = -1.0;
   H = _INFINITY;
   loopSize1 = ii - i - 1;
   loopSize2 = jj - j - 1;
   if(ii < jj) {
      N = ((2 * i)/2);
      N_loop = N;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
   } else {
      N = ((2 * j)/2);
      N_loop = 2 * jj;
      if(loopSize1 > 2) N_loop -= (loopSize1 - 2);
      if(loopSize2 > 2) N_loop -= (loopSize2 - 2);
      N_loop = (N_loop/2) - 1;
   }
#ifdef DEBUG
   if (ii <= i){
      fputs("Error in calc_bulge_internal(): ii is not greater than i\n", stderr);
   }
   if (jj <= j)
     fputs("Error in calc_bulge_internal(): jj is not greater than j\n", stderr);
#endif

#ifdef DEBUG
   if (loopSize1 + loopSize2 > maxLoop) {
      fputs("Error: calc_bulge_internal() called with loopSize1 + loopSize2 > maxLoop\n", stderr);
      free(SH);
      return;
   }
#endif
#ifdef DEBUG
   if (loopSize1 == 0 && loopSize2 == 0) {
      fputs("Error: calc_bulge_internal() called with nonsense\n", stderr);
      free(SH);
      return;
   }
#endif
   loopSize = loopSize1 + loopSize2-1;
   if((loopSize1 == 0 && loopSize2 > 0) || (loopSize2 == 0 && loopSize1 > 0)) { /* only bulges have to be considered */
      if(loopSize2 == 1 || loopSize1 == 1) { /* bulge loop of size one is treated differently
                                              the intervening nn-pair must be added */

         if((loopSize2 == 1 && loopSize1 == 0) || (loopSize2 == 0 && loopSize1 == 1)) {
            H = bulgeLoopEnthalpies[loopSize] +
              stackEnthalpies[numSeq1[i]][numSeq1[ii+off1]][numSeq2[j]][numSeq2[jj+off2]];
            S = bulgeLoopEntropies[loopSize] +
              stackEntropies[numSeq1[i]][numSeq1[ii+off1]][numSeq2[j]][numSeq2[jj+off2]];
         }
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
      } else { /* we have _not_ implemented Jacobson-Stockaymayer equation; the maximum bulgeloop size is 30 */

         H = bulgeLoopEnthalpies[loopSize] + atPenaltyH(numSeq1[i], numSeq2[j]) + atPenaltyH(numSeq1[ii+off1], numSeq2[jj+off2]);

         S = bulgeLoopEntropies[loopSize] + atPenaltyS(numSeq1[i], numSeq2[j]) + atPenaltyS(numSeq1[ii+off1], numSeq2[jj+off2]);
         if(!isFinite(H)) {
            H = _INFINITY;
            S = -1.0;
         }
        
        EntropyEnthalpy[0] = S;
        EntropyEnthalpy[1] = H;
         
      }
   } else if (loopSize1 == 1 && loopSize2 == 1) {
      S = stackint2Entropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Entropies[numSeq2[jj+off2]][numSeq2[jj+off2-1]][numSeq1[ii+off1]][numSeq1[ii+off1-1]];

      H = stackint2Enthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        stackint2Enthalpies[numSeq2[jj+off2]][numSeq2[jj+off2-1]][numSeq1[ii+off1]][numSeq1[ii+off1-1]];
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
      EntropyEnthalpy[0] = S;
      EntropyEnthalpy[1] = H;
      free(SH);
      return;
   } else { /* only internal loops */
      H = interiorLoopEnthalpies[loopSize] + tstackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEnthalpies[numSeq2[jj+off2]][numSeq2[jj+off2-1]][numSeq1[ii+off1]][numSeq1[ii+off1-1]]
        + (ILAH * abs(loopSize1 - loopSize2));

      S = interiorLoopEntropies[loopSize] + tstackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j+1]] +
        tstackEntropies[numSeq2[jj+off2]][numSeq2[jj+off2-1]][numSeq1[ii+off1]][numSeq1[ii+off1-1]] + (ILAS * abs(loopSize1 - loopSize2));
      if(!isFinite(H)) {
         H = _INFINITY;
         S = -1.0;
      }
      EntropyEnthalpy[0] = S;
      EntropyEnthalpy[1] = H;
   }
   free(SH);
   return;
}

static double 
Ss(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return -1.0;
      if (i == len1 || j == len2 + 1)
        return -1.0;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      return stackEntropies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
   } else {
      return stackEntropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}


static double 
Hs(int i, int j, int k)
{
   if(k==2) {
      if (i >= j)
        return _INFINITY;
      if (i == len1 || j == len2 + 1)
        return _INFINITY;

      if (i > len1)
        i -= len1;
      if (j > len2)
        j -= len2;
      if(isFinite(stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]])) {
         return stackEnthalpies[numSeq1[i]][numSeq1[i+1]][numSeq2[j]][numSeq2[j-1]];
      } else {
         return _INFINITY;
      }
   } else {
      return stackEnthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   }
}




/* Set default args */
void 
set_thal_default_args(thal_args *a)
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.8; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}

/* Set default args for oligo */
void
set_thal_oligo_default_args(thal_args *a)    
{
   memset(a, 0, sizeof(*a));
   a->type = thal_any; /* thal_alignment_type THAL_ANY */
   a->maxLoop = MAX_LOOP;
   a->mv = 50; /* mM */
   a->dv = 0.0; /* mM */
   a->dntp = 0.0; /* mM */
   a->dna_conc = 50; /* nM */
   a->temp = TEMP_KELVIN; /* Kelvin */
   a->dimer = 1; /* by default dimer structure is calculated */
}


static unsigned char 
/* A-0, C-1, G-2, T-3 */
str2int(char c)
{
   switch (c) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
	default:
	  return 4;
   }
   return 4;
}
/* memory stuff */

static double* 
safe_recalloc(double* ptr, int m, int n, thal_results* o)
{
   return (double*) safe_realloc(ptr, m * n * sizeof(double), o);
}

static void* 
safe_calloc(size_t m, size_t n, thal_results *o)
{
   void* ptr;
   if (!(ptr = calloc(m, n))) {
#ifdef DEBUG
      fputs("Error in calloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}

static void* 
safe_malloc(size_t n, thal_results *o)
{
   void* ptr;
   if (!(ptr = malloc(n))) {
#ifdef DEBUG
      fputs("Error in malloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}

static void* 
safe_realloc(void* ptr, size_t n, thal_results *o)
{
   ptr = realloc(ptr, n);
   if (ptr == NULL) {
#ifdef DEBUG
      fputs("Error in realloc()\n", stderr);
#endif
      THAL_OOM_ERROR;
   }
   return ptr;
}


static int 
max5(double a, double b, double c, double d, double e)
{
   if(a > b && a > c && a > d && a > e) return 1;
   else if(b > c && b > d && b > e) return 2;
   else if(c > d && c > e) return 3;
   else if(d > e) return 4;
   else return 5;
}

static void 
reverse(unsigned char *s)
{
   int i,j;
   char c;
   for (i = 0, j = length_unsig_char(s)-1; i < j; i++, j--) {
      c = s[i];
      s[i] = s[j];
      s[j] = c;
   }
}

static int 
length_unsig_char(const unsigned char * str)
{
   int i = 0;
   while(*(str++)) {
      i++;
      if(i == INT_MAX)
        return -1;
   }
   return i;
}


#define INIT_BUF_SIZE 1024

static char* 
readParamFile(const char* dirname, const char* fname, thal_results* o)
{
  FILE* file;
  char* ret = NULL;
  char* paramdir = NULL;
  paramdir = (char*) safe_malloc(strlen(dirname) + strlen(fname) + 2, o);
  strcpy(paramdir, dirname);
#ifdef OS_WIN
  if (paramdir[strlen(paramdir) - 1] != '\\') {
    strcat(paramdir, "\\\0");
  }
#else
  if (paramdir[strlen(paramdir) - 1] != '/') {
    strcat(paramdir, "/\0");
  }
#endif
  strcat(paramdir, fname);
  if (!(file = fopen(paramdir, "r"))) {
    sprintf(o->msg, "Unable to open file %s", paramdir);
    if (paramdir != NULL) {
      free(paramdir);
      paramdir = NULL;
    }
    longjmp(_jmp_buf, 1);
    return NULL;
  }
  if (paramdir != NULL) {
    free(paramdir);
    paramdir = NULL;
  }
  char c;
  int i = 0;
  size_t ssz = INIT_BUF_SIZE;
  size_t remaining_size;
  remaining_size = ssz;
  ret = (char*) safe_malloc(ssz, o);
  while (1) {
    if (feof(file)) {
      ret[i] = '\0';
      fclose(file);
      return ret;
    }
    c = fgetc(file);
    remaining_size -= sizeof(char);
    if (remaining_size <= 0) {
      if (ssz >= INT_MAX / 2) {
        strcpy(o->msg, "Out of memory");
        free(ret);
        longjmp(_jmp_buf, 1);
        return NULL;
      } else {
        ssz += INIT_BUF_SIZE;
        remaining_size += INIT_BUF_SIZE;
      }
      ret = (char *) safe_realloc(ret, ssz, o);
    }
    ret[i] = c;
    i++;
  }
}

int
thal_load_parameters(const char *path, thal_parameters *a, thal_results* o)
{
  thal_free_parameters(a);
  if (setjmp(_jmp_buf) != 0) {
    printf("longjump\n");
    return -1;
  }
  a->dangle_dh = readParamFile(path, "dangle.dh", o);
  a->dangle_ds = readParamFile(path, "dangle.ds", o);
  a->loops_dh = readParamFile(path, "loops.dh", o);
  a->loops_ds = readParamFile(path, "loops.ds", o);
  a->stack_dh = readParamFile(path, "stack.dh", o);
  a->stack_ds = readParamFile(path, "stack.ds", o);
  a->stackmm_dh = readParamFile(path, "stackmm.dh", o);
  a->stackmm_ds = readParamFile(path, "stackmm.ds", o);
  a->tetraloop_dh = readParamFile(path, "tetraloop.dh", o);
  a->tetraloop_ds = readParamFile(path, "tetraloop.ds", o);
  a->triloop_dh = readParamFile(path, "triloop.dh", o);
  a->triloop_ds = readParamFile(path, "triloop.ds", o);
  a->tstack_tm_inf_ds = readParamFile(path, "tstack_tm_inf.ds", o);
  a->tstack_dh = readParamFile(path, "tstack.dh", o);
  a->tstack2_dh = readParamFile(path, "tstack2.dh", o);
  a->tstack2_ds = readParamFile(path, "tstack2.ds", o);
  return 0;
}

static double 
saltCorrectS (double mv, double dv, double dntp)
{
   if(dv<=0) dntp=dv;
   return 0.368*((log((mv+120*(sqrt(fmax(0.0, dv-dntp))))/1000)));
}

static char*
th_read_str_line(char **str, thal_results* o)
{
  if (*str == NULL) {
    return NULL;
  }
  char *ptr = *str;
  char *ini = *str;
  while(1) {
    if ((*ptr == '\n') || (*ptr == '\0')) {
      char *ret = NULL;
      if (!(ret = malloc(sizeof(char) * (ptr - ini + 1)))) {
#ifdef DEBUG
        fputs("Error in malloc()\n", stderr);
#endif
        THAL_OOM_ERROR;
      }
      /* copy line */
      strncpy(ret, ini, (ptr - ini + 1));
      ret[ptr - ini] = '\0';

      if (*ptr == '\0') { /* End of String */
        *str = NULL;
      } else {
        ptr++;
        if (*ptr == '\0') { /* End of String */
          *str = NULL;
        } else {
          *str = ptr;
        }
      }
      if (ptr == ini) {
        if (ret != NULL) {
          free(ret);
        }
        return NULL;
      } else {  
        return ret;
      }
    }
    ptr++;
  }
}

/* These functions are needed as "inf" cannot be read on Windows directly */
static double 
readDouble(char **str, thal_results* o)
{
  double result;
  char *line = th_read_str_line(str, o);
  /* skip any spaces at beginning of the line */
  while (isspace(*line)) line++;
  if (!strncmp(line, "inf", 3)) {
    free(line);
    return _INFINITY;
  }
  sscanf(line, "%lf", &result);
  if (line != NULL) {
    free(line);
  }
  return result;
}

/* Reads a line containing 4 doubles, which can be specified as "inf". */
static void
readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o)
{
  char *line = th_read_str_line(str, o);
  char *p = line, *q;
  /* skip first number on the line */
  while (isspace(*p)) p++;
  while (isdigit(*p)) p++;
  while (isspace(*p)) p++;
  /* read second number */
  q = p;
  while (!isspace(*q)) q++;
  *q = '\0'; q++;
  if (!strcmp(p, "inf")) *v1 = _INFINITY;
  else sscanf(p, "%lf", v1);
  while (isspace(*q)) q++;
  /* read third number */
  p = q;
  while (!isspace(*p)) p++;
  *p = '\0'; p++;
  if (!strcmp(q, "inf")) *v2 = _INFINITY;
  else sscanf(q, "%lf", v2);
  while (isspace(*p)) p++;
  /* read last number */
  q = p;
  while (!isspace(*q) && (*q != '\0')) q++;
  *q = '\0';
  if (!strcmp(p, "inf")) *v3 = _INFINITY;
  else sscanf(p, "%lf", v3);
  if (line != NULL) {
    free(line);
  }
}

/* Reads a line containing a short string and a double, used for reading a triloop or tetraloop. */
static int
readTLoop(char **str, char *s, double *v, int triloop, thal_results *o)
{
  char *line = th_read_str_line(str, o);
  if (!line) return -1;
  char *p = line, *q;
  /* skip first spaces */
  while (isspace(*p)) p++;
  /* read the string */
  q = p;
  while (isalpha(*q)) q++;
  *q = '\0'; q++;
  if (triloop) {
    strncpy(s, p, 5);   /*triloop string has 5 characters*/
    s[5] = '\0';
  } else {
    strncpy(s, p, 6);   /*tetraloop string has 6 characters*/
    s[6] = '\0';
  }
  /* skip all spaces */
  while (isspace(*q)) q++;
  p = q;
  while (!isspace(*p) && (*p != '\0')) p++;
  *p = '\0';
  if (!strcmp(q, "inf")) *v = _INFINITY;
  else sscanf(q, "%lg", v);
  if (line != NULL) {
    free(line);
  }
  return 0;
} 

static void 
getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stack_ds;
   char *pt_dh = tp->stack_dh;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackEntropies[i][ii][j][jj] = -1.0;
                  stackEnthalpies[i][ii][j][jj] = _INFINITY;
               } else {
                  stackEntropies[i][ii][j][jj] = readDouble(&pt_ds, o);
                  stackEnthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
                  if (!isFinite(stackEntropies[i][ii][j][jj]) || !isFinite(stackEnthalpies[i][ii][j][jj])) {
                     stackEntropies[i][ii][j][jj] = -1.0;
                     stackEnthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

static void 
getStackint2(double stackint2Entropies[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, ii, jj;
   char *pt_ds = tp->stackmm_ds;
   char *pt_dh = tp->stackmm_dh;
   int N;
   for (i = 0; i < 5; ++i) {
      for (ii = 0; ii < 5; ++ii) {
         for (j = 0; j < 5; ++j) {
            for (jj = 0; jj < 5; ++jj) {
               if (i == 4 || j == 4 || ii == 4 || jj == 4) {
                  stackint2Entropies[i][ii][j][jj] = -1.0;
                  stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
               } else {
				   N++;
                  stackint2Entropies[i][ii][j][jj] = readDouble(&pt_ds, o);
                  stackint2Enthalpies[i][ii][j][jj] = readDouble(&pt_dh, o);
//				   printf("N=%d, i=%d, ii=%d, j=%d, jj=%d, result=%f\n", N, i, ii, j, jj, stackint2Enthalpies[i][ii][j][jj]);
                  if (!isFinite(stackint2Entropies[i][ii][j][jj]) || !isFinite(stackint2Enthalpies[i][ii][j][jj])) {
                     stackint2Entropies[i][ii][j][jj] = -1.0;
                     stackint2Enthalpies[i][ii][j][jj] = _INFINITY;
                  }
               }
            }
         }
      }
   }
}

/*
static void 
verifyStackTable(double stack[5][5][5][5], char* type)
{

   int i, j, ii, jj;
   for (i = 0; i < 4; ++i)
     for (j = 0; j < 4; ++j)
       for (ii = 0; ii < 4; ++ii)
         for (jj = 0; jj < 4; ++jj)
           if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
#ifdef DEBUG
             fprintf(stderr, "Warning: symmetrical stacks _are_ _not_ equal: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n",
#endif
                     BASES[i], BASES[j], BASES[ii], BASES[jj], type, stack[i][j][ii][jj], BASES[jj],
                     BASES[ii], BASES[j], BASES[i], type, stack[jj][ii][j][i]);
}
*/

static void 
getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
          double dangleEnthalpies5[5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i, j, k;
   char *pt_ds = tp->dangle_ds;
   char *pt_dh = tp->dangle_dh;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies3[i][k][j] = -1.0;
             dangleEnthalpies3[i][k][j] = _INFINITY;
          } else {
             dangleEntropies3[i][k][j] = readDouble(&pt_ds, o);
             dangleEnthalpies3[i][k][j] = readDouble(&pt_dh, o);
             if(!isFinite(dangleEntropies3[i][k][j]) || !isFinite(dangleEnthalpies3[i][k][j])) {
                dangleEntropies3[i][k][j] = -1.0;
                dangleEnthalpies3[i][k][j] = _INFINITY;             
             }
          }
       }

   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       for (k = 0; k < 5; ++k) {
          if (i == 4 || j == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else if (k == 4) {
             dangleEntropies5[i][j][k] = -1.0;
             dangleEnthalpies5[i][j][k] = _INFINITY;
          } else {
             dangleEntropies5[i][j][k] = readDouble(&pt_ds, o);
             dangleEnthalpies5[i][j][k] = readDouble(&pt_dh, o);
             if(!isFinite(dangleEntropies5[i][j][k]) || !isFinite(dangleEnthalpies5[i][j][k])) {
                dangleEntropies5[i][j][k] = -1.0;
                dangleEnthalpies5[i][j][k] = _INFINITY;
             }
          }
       }
}

static void 
getLoop(double hairpinLoopEntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropies[30],
        double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], 
        const thal_parameters *tp, thal_results* o)
{
   int k;
   char *pt_ds = tp->loops_ds;
   char *pt_dh = tp->loops_dh;
   for (k = 0; k < 30; ++k) {
      readLoop(&pt_ds, &interiorLoopEntropies[k], &bulgeLoopEntropies[k], &hairpinLoopEntropies[k], o);
      readLoop(&pt_dh, &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k], o);
   }
}

static void 
getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{
   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack_tm_inf_ds;
   char *pt_dh = tp->tstack_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4) {
              tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              tstackEntropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstackEntropies[i1][i2][j1][j2] = 0.00000000001;
              tstackEnthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstackEntropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
              tstackEnthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
              if (!isFinite(tstackEntropies[i1][i2][j1][j2]) || !isFinite(tstackEnthalpies[i1][i2][j1][j2])) {
                 tstackEntropies[i1][i2][j1][j2] = -1.0;
                 tstackEnthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void 
getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results* o)
{

   int i1, j1, i2, j2;
   char *pt_ds = tp->tstack2_ds;
   char *pt_dh = tp->tstack2_dh;
   for (i1 = 0; i1 < 5; ++i1)
     for (i2 = 0; i2 < 5; ++i2)
       for (j1 = 0; j1 < 5; ++j1)
         for (j2 = 0; j2 < 5; ++j2)
           if (i1 == 4 || j1 == 4)  {
              tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              tstack2Entropies[i1][i2][j1][j2] = -1.0;
           } else if (i2 == 4 || j2 == 4) {
              tstack2Entropies[i1][i2][j1][j2] = 0.00000000001;
              tstack2Enthalpies[i1][i2][j1][j2] = 0.0;
           } else {
              tstack2Entropies[i1][i2][j1][j2] = readDouble(&pt_ds, o);
              tstack2Enthalpies[i1][i2][j1][j2] = readDouble(&pt_dh, o);
              if (!isFinite(tstack2Entropies[i1][i2][j1][j2]) || !isFinite(tstack2Enthalpies[i1][i2][j1][j2])) {
                 tstack2Entropies[i1][i2][j1][j2] = -1.0;
                 tstack2Enthalpies[i1][i2][j1][j2] = _INFINITY;
              }
           }
}

static void 
tableStartATS(double atp_value, double atpS[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpS[i][j] = 0.00000000001;
   atpS[0][3] = atpS[3][0] = atp_value;
}


static void 
tableStartATH(double atp_value, double atpH[5][5])
{

   int i, j;
   for (i = 0; i < 5; ++i)
     for (j = 0; j < 5; ++j)
       atpH[i][j] = 0.0;

   atpH[0][3] = atpH[3][0] = atp_value;
}



/* Initialize the thermodynamic values (parameters) */
int  thal_set_null_parameters(thal_parameters *a) {
  a->dangle_dh = NULL;
  a->dangle_ds = NULL;
  a->loops_dh = NULL;
  a->loops_ds = NULL;
  a->stack_dh = NULL;
  a->stack_ds = NULL;
  a->stackmm_dh = NULL;
  a->stackmm_ds = NULL;
  a->tetraloop_dh = NULL;
  a->tetraloop_ds = NULL;
  a->triloop_dh = NULL;
  a->triloop_ds = NULL;
  a->tstack_tm_inf_ds = NULL;
  a->tstack_dh = NULL;
  a->tstack2_dh = NULL;
  a->tstack2_ds = NULL;
  return 0;
}

/* Free the thermodynamic values (parameters) */
int  thal_free_parameters(thal_parameters *a) {
  if (NULL != a->dangle_dh) {
    free(a->dangle_dh);
    a->dangle_dh = NULL;
  }
  if (NULL != a->dangle_ds) {
    free(a->dangle_ds);
    a->dangle_ds = NULL;
  }
  if (NULL != a->loops_dh) {
    free(a->loops_dh);
    a->loops_dh = NULL;
  }
  if (NULL != a->loops_ds) {
    free(a->loops_ds);
    a->loops_ds = NULL;
  }
  if (NULL != a->stack_dh) {
    free(a->stack_dh);
    a->stack_dh = NULL;
  }
  if (NULL != a->stack_ds) {
    free(a->stack_ds);
    a->stack_ds = NULL;
  }
  if (NULL != a->stackmm_dh) {
    free(a->stackmm_dh);
    a->stackmm_dh = NULL;
  }
  if (NULL != a->stackmm_ds) {
    free(a->stackmm_ds);
    a->stackmm_ds = NULL;
  }
  if (NULL != a->tetraloop_dh) {
    free(a->tetraloop_dh);
    a->tetraloop_dh = NULL;
  }
  if (NULL != a->tetraloop_ds) {
    free(a->tetraloop_ds);
    a->tetraloop_ds = NULL;
  }
  if (NULL != a->triloop_dh) {
    free(a->triloop_dh);
    a->triloop_dh = NULL;
  }
  if (NULL != a->triloop_ds) {
    free(a->triloop_ds);
    a->triloop_ds = NULL;
  }
  if (NULL != a->tstack_tm_inf_ds) {
    free(a->tstack_tm_inf_ds);
    a->tstack_tm_inf_ds = NULL;
  }
  if (NULL != a->tstack_dh) {
    free(a->tstack_dh);
    a->tstack_dh = NULL;
  }
  if (NULL != a->tstack2_dh) {
    free(a->tstack2_dh);
    a->tstack2_dh = NULL;
  }
  if (NULL != a->tstack2_ds) {
    free(a->tstack2_ds);
    a->tstack2_ds = NULL;
  }
  return 0;
}




/* Read the thermodynamic values (parameters) from the parameter files
   in the directory specified by 'path'.  Return 0 on success and -1
   on error. The thermodynamic values are stored in multiple static
   variables. */
int 
get_thermodynamic_values(const thal_parameters *tp, thal_results *o)
{
  if (setjmp(_jmp_buf) != 0) {
     return -1;
  }
  getStack(stackEntropies, stackEnthalpies, tp, o);
  /* verifyStackTable(stackEntropies, "entropy");
     verifyStackTable(stackEnthalpies, "enthalpy"); */ /* this is for code debugging */
  getStackint2(stackint2Entropies, stackint2Enthalpies, tp, o);
  getDangle(dangleEntropies3, dangleEnthalpies3, dangleEntropies5, dangleEnthalpies5, tp, o);
  getLoop(hairpinLoopEntropies, interiorLoopEntropies, bulgeLoopEntropies, hairpinLoopEnthalpies,
          interiorLoopEnthalpies, bulgeLoopEnthalpies, tp, o);
  getTstack(tstackEntropies, tstackEnthalpies, tp, o);
  getTstack2(tstack2Entropies, tstack2Enthalpies, tp, o);
  /* getting the AT-penalties */
  tableStartATS(AT_S, atpS);
  tableStartATH(AT_H, atpH);

  return 0;
}

void 
destroy_thal_structures()
{
  if (triloopEntropies != NULL) {
    free(triloopEntropies);
    triloopEntropies = NULL;
  }
  if (triloopEnthalpies != NULL) {
    free(triloopEnthalpies);
    triloopEnthalpies = NULL;
  }
  if (tetraloopEntropies != NULL) {
    free(tetraloopEntropies);
    tetraloopEntropies = NULL;
  }
  if (tetraloopEnthalpies != NULL) {
    free(tetraloopEnthalpies);
    tetraloopEnthalpies = NULL;
  }
}

static void 
RSH(int i, int j, double* EntropyEnthalpy)
{
   double G1, G2;
   double S1, S2;
   double H1, H2;
   double T1, T2;
   S1 = S2 = -1.0;
   H1 = H2 = _INFINITY;
   T1 = T2 = -_INFINITY;
   if (bpIndx(numSeq1[i], numSeq2[j]) == 0) {
      EntropyEnthalpy[0] = -1.0;
      EntropyEnthalpy[1] = _INFINITY;
      return;
   }
   S1 = atPenaltyS(numSeq1[i], numSeq2[j]) + tstack2Entropies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   H1 = atPenaltyH(numSeq1[i], numSeq2[j]) + tstack2Enthalpies[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]][numSeq2[j + 1]];
   G1 = H1 - TEMP_KELVIN*S1;
   if(!isFinite(H1) || G1>0) {
      H1 = _INFINITY;
      S1 = -1.0;
      G1 = 1.0;
   }
   
   if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]]) && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]] +
        dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
    G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies3[numSeq1[i]][numSeq1[i + 1]][numSeq2[j]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2 >0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if(G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }

   else if(bpIndx(numSeq1[i+1], numSeq2[j+1]) == 0 && isFinite(dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]])) {
      S2 = atPenaltyS(numSeq1[i], numSeq2[j]) + dangleEntropies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      H2 = atPenaltyH(numSeq1[i], numSeq2[j]) + dangleEnthalpies5[numSeq1[i]][numSeq2[j]][numSeq2[j + 1]];
      G2 = H2 - TEMP_KELVIN*S2;
      if(!isFinite(H2) || G2>0) {
         H2 = _INFINITY;
         S2 = -1.0;
     G2 = 1.0;
      }
      T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
      if(isFinite(H1) && G1<0) {
         T1 = (H1 + dplx_init_H) / (S1 + dplx_init_S + RC);
         if(T1 < T2 && G2<0) {
            S1 = S2;
            H1 = H2;
            T1 = T2;
         }
      } else if (G2<0){
         S1 = S2;
         H1 = H2;
         T1 = T2;
      }
   }
   S2 = atPenaltyS(numSeq1[i], numSeq2[j]);
   H2 = atPenaltyH(numSeq1[i], numSeq2[j]);
   T2 = (H2 + dplx_init_H) / (S2 + dplx_init_S + RC);
   G1 = H1 -TEMP_KELVIN*S1;
   G2 =  H2 -TEMP_KELVIN*S2;
   if(isFinite(H1)) {
      if(T1 < T2) {
         EntropyEnthalpy[0] = S2;
         EntropyEnthalpy[1] = H2;
      } else {
         EntropyEnthalpy[0] = S1;
         EntropyEnthalpy[1] = H1;
      }
   } else {
      EntropyEnthalpy[0] = S2;
      EntropyEnthalpy[1] = H2;
   }
   return;
}


