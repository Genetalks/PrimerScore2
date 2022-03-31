#include <cassert>
#include "tm_calculator.h"
#include "alignment_TM.h"


pt::tm_calculator::tm_calculator(const pt::options &opt){
	o.sec_struct = NULL;
	set_thal_default_args(&a);
	a.mv = opt.mv;
	a.dv = opt.dv;
	a.dna_conc = opt.dna;
	a.temp = opt.temp + ABSOLUTE_ZERO;
	a.dntp = opt.dntp;
	a.maxLoop = MAX_LOOP;
	thal_set_null_parameters(&thermodynamic_parameters);

	thal_load_parameters(opt.path.c_str(), &thermodynamic_parameters, &o);
	get_thermodynamic_values(&thermodynamic_parameters, &o);

	/*** Calc part of the salt correction ***/
	saltCorrection = saltCorrectS(a.mv, a.dv, a.dntp); /* salt correction for entropy, must be multiplied with N, which is
														the total number of phosphates in the duplex divided by 2; 8bp dplx N=7 */
}

pt::tm_calculator::~tm_calculator(){
   destroy_thal_structures();
   thal_free_parameters(&thermodynamic_parameters);
}


void pt::tm_calculator::alignment_tm(const std::string &align0, const std::string &seq1, const std::string &seq2, double *tm, double *dG, double *dS, double *dH) {
	assert(align0.size() == seq1.size());
	assert(align0.size() == seq2.size());

	int i, j;
	int len = align0.size();
	double H, S, G, t;
	double *SH;
	int end1, end2;
	dplx_init_H = 200;
	dplx_init_S = -5.7;
	RC = R * log(a.dna_conc / 4000000000.0);

	unsigned char *numSeq1, *numSeq2, *align; /* same as seq1 and seq2 but converted to numbers */
	/* convert nucleotides to numbers */
	numSeq1 = (unsigned char *)safe_malloc(len + 2, &this->o);
	numSeq2 = (unsigned char *)safe_malloc(len + 2, &this->o);
	align = (unsigned char *)safe_malloc(len + 2, &this->o);

	// reverse(seq2); /* REVERSE oligo2, so it goes to dpt 3'->5' direction */
	for (i = 1; i <= len; ++i)
	{
		if (align0[i - 1] == '#')
		{
			numSeq1[i] = 0; // tm of 0+1(numSeq1=0, numSeq2=1) is 1 less than tm of 0+0;
			numSeq2[i] = 1; // tm of 1+1(numSeq1=1, numSeq2=1) is 2-3 less than tm of 0+0;
		}
		else
		{
			/* A-0, C-1, G-2, T-3 */
			numSeq1[i] = str2int(seq1[i - 1]);
			// numSeq2[i] = str2int(seq2[i - 1]);
			numSeq2[i] = str2int(seq2[len - i]);
		}
		align[i] = align0[i - 1];
	}
	numSeq1[0] = numSeq1[len + 1] = 4; /* mark as end */
	numSeq2[0] = numSeq2[len + 1] = 4; /* mark as N-s */
	align[0] = align[len + 1] = '\0';

	/*** caculate tm ***/
	S = dplx_init_S;
	H = dplx_init_H;
	SH = (double *)safe_malloc(2 * sizeof(double), &this->o);
	int iterminal = -1;
	int N = 0;
	for (i = 1; i <= len; i++)
	{
		if (align[i] != '|')
		{
			continue;
		}
		N++;
		if (iterminal == -1)
		{ /* Terminal AT penalty */
			S += atPenaltyS(numSeq1[i], numSeq2[i]);
			H += atPenaltyH(numSeq1[i], numSeq2[i]);
		}
		iterminal = i;

		/* left */
		if (align[i - 1] == '#')
		{ /* left dange */
			SH_Ldange(numSeq1, numSeq2, i, i, SH);
			S += SH[0];
			H += SH[1];
			//		}else if(align[i-1] == '*'){ /* left mismatch */
			//			S += tstack2Entropies[numSeq2[i]][numSeq2[i-1]][numSeq1[i]][numSeq1[i-1]];
			//			H += tstack2Enthalpies[numSeq2[i]][numSeq2[i-1]][numSeq1[i]][numSeq1[i-1]];
		}
		else if (align[i - 1] == '|')
		{ /* left align */
			S += stackEntropies[numSeq1[i - 1]][numSeq1[i]][numSeq2[i - 1]][numSeq2[i]];
			H += stackEnthalpies[numSeq1[i - 1]][numSeq1[i]][numSeq2[i - 1]][numSeq2[i]];
		}

		if (align[i + 1] == '#')
		{ /* right dange */
			SH_Rdange(numSeq1, numSeq2, i, i, SH);
			S += SH[0];
			H += SH[1];
			//		}else if(align[i+1]== '^' || align[i+1] =='-' || (align[i+1]=='*' && i+2<len && align[i+2]=='*')){ /* loop */
		}
		else if (align[i + 1] == '^' || align[i + 1] == '-' || align[i + 1] == '*')
		{ /* loop */
			end1 = end2 = i;
			int off1 = 0;
			int off2 = 0;
			for (j = i + 1; j < len + 1; j++)
			{
				if (align[j] == '^')
				{
					end1++;
					off2++;
				}
				else if (align[j] == '-')
				{
					off1++;
					end2++;
				}
				else if (align[j] == '*')
				{
					end1++;
					end2++;
				}
				else if (align[j] == '|')
				{
					end1++;
					end2++;
					break;
				}
			}
			calc_bulge_internal(numSeq1, numSeq2, i, i, end1, end2, SH, 0, a.maxLoop, off1, off2);
			i = j - 1;
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
	G = (H) - ((ABSOLUTE_ZERO + 37) * (S + (N * saltCorrection)));
	S = S + (N * saltCorrection);
	if(H>0){
		t = -100;
	}

	// printf("dS=%f, dH=%f, dG=%f, t=%f\n", S, H, G, t);
	*dS = S;
	*dH = H;
	*dG = G;
	*tm = t;

	free(numSeq1);
	free(numSeq2);
	free(align);
	free(SH);
}
