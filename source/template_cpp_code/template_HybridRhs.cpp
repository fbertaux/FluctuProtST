/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "HybridRhs.hpp"

HybridRhs::HybridRhs ( ModelParameters* modelParameters ) : mf_ModelParameters (modelParameters) {}


void
HybridRhs::setMrnaTables(MatDoub *mrtab, VecDoub *ttab, int neventss)
{
    mrnatable=mrtab;
    timetable=ttab;
    tindex=0;
    nevents=neventss;
    mf_computedReactionRates = VecDoub ( mf_ModelParameters->mf_NumReacs , 0. ) ;
}


void
HybridRhs::findGoodTindex(const Doub x)
{
    if (tindex<0) { cout << "NEGATIVE TINDEX !!!" << endl; exit(234); }
    if (tindex>nevents-1) { cout << "TOO BIG TINDEX !!!" << endl; exit(235); }
    while ( ! ( ( x >= (*timetable)[tindex] ) && ( x <= (*timetable)[tindex+1] ) ) )
    {
        if ( x > (*timetable)[tindex+1] )
        {
            if (tindex <= nevents - 2) { tindex++; }
//            else { tindex = nevents - 1 ; return ; } // "unsafe" version
            else { cout << "\t\t\t\tmrna time in the future asked? I quit." << endl; exit(236); }
        }
        else
        {
            if (tindex>0) { tindex--; }
            else { cout << "\t\t\t\tmrna time in the past asked? I quit." << endl; exit(237); }
        }
    }
}


void
HybridRhs::operator () (const Doub x, VecDoub_I &y, VecDoub_O &dydx)
{
    findGoodTindex (x) ;

 placeholder_hybrid_rhs}


