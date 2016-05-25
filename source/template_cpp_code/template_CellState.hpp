/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "libs/nr3.h"

#include "ModelParameters.hpp"

struct CellState
{
	ModelParameters* mf_ModelParameters ;
	VecDoub mf_AllProts ;
	VecDoub mf_GeneMrnas ;

	CellState ( ModelParameters* modelParameters ) ;

	// methods for name access to species
placeholder_name_access};

