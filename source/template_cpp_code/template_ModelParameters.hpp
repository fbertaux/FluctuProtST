/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#ifndef MODEL_PARAMETERS
#define MODEL_PARAMETERS

#include "libs/nr3.h"

struct ModelParameters
{
	Int mf_NumGenes ;
	Int mf_NumAllProteinSpecies ;
	Int mf_NumModifiedProteins ;
	Int mf_NumReacs ;
	VecDoub mf_kons ;
	VecDoub mf_koffs ;
	VecDoub mf_ksms ;
	VecDoub mf_rms ;
	VecDoub mf_ksps ;
	VecDoub mf_rps ;
	VecDoub mf_kreacs ;
	VecDoub mf_degrates ;

	ModelParameters () ;
};

#endif
