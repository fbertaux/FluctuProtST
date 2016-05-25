/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "ModelParameters.hpp"

struct HybridRhs
{
    HybridRhs ( ModelParameters* modelParameters ) ;
	ModelParameters* mf_ModelParameters ;


	// fields 
    int tindex;
    MatDoub *mrnatable;
    VecDoub *timetable;
    int nevents;
    VecDoub mf_computedReactionRates ;

    // methods
    void setMrnaTables (MatDoub *mrtab, VecDoub* ttab, int neventss);
    void findGoodTindex (const Doub x);
    void operator() (const Doub x, VecDoub_I &y, VecDoub_O &dydx);

};

