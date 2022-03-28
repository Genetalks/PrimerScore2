#ifndef _ALIGNMENT_TM_H
#define _ALIGNMENT_TM_H

#include <limits.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "thal.h"

#ifdef __cplusplus
extern "C"
{
#endif

#define CHECK_ERROR(COND, MSG) \
  if (COND)                    \
  {                            \
    strcpy(o->msg, MSG);       \
    errno = 0;                 \
    longjmp(_jmp_buf, 1);      \
  }
#define THAL_OOM_ERROR               \
  {                                  \
    strcpy(o->msg, "Out of memory"); \
    errno = ENOMEM;                  \
    longjmp(_jmp_buf, 1);            \
  }
#define THAL_IO_ERROR(f)                          \
  {                                               \
    sprintf(o->msg, "Unable to open file %s", f); \
    longjmp(_jmp_buf, 1);                         \
  }

#define bpIndx(a, b) BPI[a][b] /* for traceing matrix BPI */
#define atPenaltyS(a, b) atpS[a][b]
#define atPenaltyH(a, b) atpH[a][b]

#define SMALL_NON_ZERO 0.000001
#define DBL_EQ(X, Y) (((X) - (Y)) < (SMALL_NON_ZERO) ? (1) : (2)) /* 1 when numbers are equal */

#ifdef INTEGER
#define isFinite(x) (x < _INFINITY / 2)
#else
#define isFinite(x) isfinite(x)
#endif

#define isPositive(x) ((x) > 0 ? (1) : (0))

  /*** BEGIN CONSTANTS ***/
  extern const double _INFINITY;

  extern const double R;                /* cal/Kmol */
  extern const double ILAS;             /* Internal Loop Entropy ASymmetry correction -0.3kcal/mol*/
  extern const double ILAH;             /* Internal Loop EntHalpy Asymmetry correction */
  extern const double AT_H;             /* AT penalty */
  extern const double AT_S;             /* AT penalty */
  extern const double MinEntropyCutoff; /* to filter out non-existing entropies */
  extern const double MinEntropy;       /* initiation */
  extern const double G2;               /* structures w higher G are considered to be unstabile */
  extern const double ABSOLUTE_ZERO;
  extern const double TEMP_KELVIN;

  // const int MAX_LOOP = 30; /* the maximum size of loop that can be calculated; for larger loops formula must be implemented */
  // const int MIN_LOOP = 0;

  // static const char BASES[5] = {'A', 'C', 'G', 'T', 'N'}; /* bases to be considered - N is every symbol that is not A, G, C,$
  //                                                   */
  // static const char BASE_PAIRS[4][4] = {"A-T", "C-G", "G-C", "T-A" }; /* allowed basepairs */
  /* matrix for allowed; bp 0 - no bp, watson crick bp - 1 */
  static const int BPI[6][6] = {
      {0, 0, 0, 1, 0, 0}, /* A, C, G, T, N, end; */
      {0, 0, 1, 0, 0, 0},
      {0, 1, 0, 0, 0, 0},
      {1, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0}, /* A, C, G, T, N, end; */
      {0, 0, 0, 0, 0, 0}};

  /*** END OF CONSTANTS ***/

  /*** BEGIN STRUCTs ***/

  struct triloop
  {
    char loop[5];
    double value;
  };

  struct tetraloop
  {
    char loop[6];
    double value;
  };

  struct tracer /* structure for tracebacku - unimolecular str */
  {
    int i;
    int j;
    int mtrx; /* [0 1] EntropyDPT/EnthalpyDPT*/
    struct tracer *next;
  };

  /*** END STRUCTs ***/

  extern void SH_Ldange(unsigned char *numSeq1, unsigned char *numSeq2, int i, int j, double *EntropyEnthalpy);
  extern void SH_Rdange(unsigned char *numSeq1, unsigned char *numSeq2, int i, int j, double *EntropyEnthalpy);
  extern unsigned char str2int(char c);                                                /* converts DNA sequence to int; 0-A, 1-C, 2-G, 3-T, 4-whatever */
  extern int length_unsig_char(const unsigned char *str);                              /* returns length of unsigned char; to avoid warnings while compiling */
  extern double saltCorrectS(double mv, double dv, double dntp);                       /* part of calculating salt correction
                                                                                          for Tm by SantaLucia et al */
  extern char *readParamFile(const char *dirname, const char *fname, thal_results *o); /* file of thermodynamic params */

  /* get thermodynamic tables */
  extern double readDouble(char **str, thal_results *o);

  extern void readLoop(char **str, double *v1, double *v2, double *v3, thal_results *o);

  extern int readTLoop(char **str, char *s, double *v, int triloop, thal_results *o);

  extern void getStack(double stackEntropies[5][5][5][5], double stackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results *o);

  extern void getStackint2(double stackEntropiesint2[5][5][5][5], double stackint2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results *o);

  extern void getDangle(double dangleEntropies3[5][5][5], double dangleEnthalpies3[5][5][5], double dangleEntropies5[5][5][5],
                        double dangleEnthalpies5[5][5][5], const thal_parameters *tp, thal_results *o);

  extern void getTstack(double tstackEntropies[5][5][5][5], double tstackEnthalpies[5][5][5][5], const thal_parameters *tp, thal_results *o);

  extern void getTstack2(double tstack2Entropies[5][5][5][5], double tstack2Enthalpies[5][5][5][5], const thal_parameters *tp, thal_results *o);

  extern void getTriloop(struct triloop **, struct triloop **, int *num, const thal_parameters *tp, thal_results *o);

  extern void getTetraloop(struct tetraloop **, struct tetraloop **, int *num, const thal_parameters *tp, thal_results *o);

  extern void getLoop(double hairpinLoopEnntropies[30], double interiorLoopEntropies[30], double bulgeLoopEntropiess[30],
                      double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], const thal_parameters *tp, thal_results *o);

  extern void tableStartATS(double atp_value, double atp[5][5]); /* creates table of entropy values for nucleotides
                                                                    to which AT-penlty must be applied */

  extern void tableStartATH(double atp_value, double atp[5][5]);

  extern int comp3loop(const void *, const void *); /* checks if sequnece consists of specific triloop */

  extern int comp4loop(const void *, const void *); /* checks if sequnece consists of specific tetraloop */

  /* calculates bulges and internal loops for dimer structures */
  extern void calc_bulge_internal(unsigned char *numSeq1, unsigned char *numSeq2, int ii, int jj, int i, int j, double *EntropyEnthalpy, int traceback, int maxLoop, int off1, int off2);

  extern double Hd5(int, int);     /* returns thermodynamic value (H) for 5' dangling end */
  extern double Hd3(int, int);     /* returns thermodynamic value (H) for 3' dangling end */
  extern double Sd5(int, int);     /* returns thermodynamic value (S) for 5' dangling end */
  extern double Sd3(int, int);     /* returns thermodynamic value (S) for 3' dangling end */
  extern double Ststack(int, int); /* returns entropy value for terminal stack */
  extern double Htstack(int, int); /* returns enthalpy value for terminal stack */

  /* memory stuff */
  extern void *safe_calloc(size_t, size_t, thal_results *o);
  extern void *safe_malloc(size_t, thal_results *o);
  extern void *safe_realloc(void *, size_t, thal_results *o);
  extern double *safe_recalloc(double *ptr, int m, int n, thal_results *o);

  extern double atpS[5][5];     /* AT penalty */
  extern double atpH[5][5];     /* AT penalty */
  extern double *send5, *hend5; /* calc 5'  */
  /* w/o init not constant anymore, cause for unimolecular and bimolecular foldings there are different values */
  extern double dplx_init_H;    /* initiation enthalpy; for duplex 200, for unimolecular structure 0 */
  extern double dplx_init_S;    /* initiation entropy; for duplex -5.7, for unimoleculat structure 0 */
  extern double saltCorrection; /* value calculated by saltCorrectS, includes correction for monovalent and divalent cations */
  extern double RC;             /* universal gas constant multiplied w DNA conc - for melting temperature */
  extern double SHleft;         /* var that helps to find str w highest melting temperature */

  extern double dangleEntropies3[5][5][5];       /* thermodynamic paramteres for 3' dangling ends */
  extern double dangleEnthalpies3[5][5][5];      /* ther params for 3' dangling ends */
  extern double dangleEntropies5[5][5][5];       /* ther params for 5' dangling ends */
  extern double dangleEnthalpies5[5][5][5];      /* ther params for 5' dangling ends */
  extern double stackEntropies[5][5][5][5];      /* ther params for perfect match pairs */
  extern double stackEnthalpies[5][5][5][5];     /* ther params for perfect match pairs */
  extern double stackint2Entropies[5][5][5][5];  /*ther params for perfect match and internal mm */
  extern double stackint2Enthalpies[5][5][5][5]; /* ther params for perfect match and internal mm*/
  extern double interiorLoopEntropies[30];       /* interior loop params according to length of the loop */
  extern double bulgeLoopEntropies[30];          /* bulge loop params according to length of the loop */
  extern double hairpinLoopEntropies[30];        /* hairpin loop params accordint to length of the loop */
  extern double interiorLoopEnthalpies[30];      /* same as interiorLoopEntropies but values of entropy */
  extern double bulgeLoopEnthalpies[30];         /* same as bulgeLoopEntropies but values of entropy */
  extern double hairpinLoopEnthalpies[30];       /* same as hairpinLoopEntropies but values of entropy */
  extern double tstackEntropies[5][5][5][5];     /* ther params for terminal mismatches */
  extern double tstackEnthalpies[5][5][5][5];    /* ther params for terminal mismatches */
  extern double tstack2Entropies[5][5][5][5];    /* ther params for internal terminal mismatches */
  extern double tstack2Enthalpies[5][5][5][5];   /* ther params for internal terminal mismatches */
  extern struct triloop *triloopEntropies;       /* ther penalties for given triloop seq-s */
  extern struct triloop *triloopEnthalpies;      /* ther penalties for given triloop seq-s */
  extern struct tetraloop *tetraloopEntropies;   /* ther penalties for given tetraloop seq-s */
  extern struct tetraloop *tetraloopEnthalpies;  /* ther penalties for given tetraloop seq-s */
  extern jmp_buf _jmp_buf;

#ifdef __cplusplus
}
#endif

#endif
