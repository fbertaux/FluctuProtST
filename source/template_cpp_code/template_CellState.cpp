/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "CellState.hpp"

CellState::CellState ( ModelParameters* modelParameters ) : mf_ModelParameters (modelParameters)
{
	mf_AllProts = VecDoub ( mf_ModelParameters->mf_NumAllProteinSpecies , 0. ) ;
	mf_GeneMrnas = VecDoub ( 3*mf_ModelParameters->mf_NumGenes , 0.) ;
    for ( Int g = 0 ; g < mf_ModelParameters->mf_NumGenes ; g++ )
    {
        Doub EG = mf_ModelParameters->mf_kons[g] / ( mf_ModelParameters->mf_kons[g] + mf_ModelParameters->mf_koffs[g] ) ;
        Doub EM = EG * mf_ModelParameters->mf_ksms[g] / mf_ModelParameters->mf_rms[g] ;
        mf_GeneMrnas[3*g] = 1. ; // gene on
        mf_GeneMrnas[3*g+2] = (Int) EM ; // mrna ~ mean
        Doub EP = EM * mf_ModelParameters->mf_ksps[g] / mf_ModelParameters->mf_rps[g] ;
        mf_AllProts[g] = EP ;
    }
}
