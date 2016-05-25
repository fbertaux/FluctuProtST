/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/


#include "CellState.hpp"

#include "libs/ran.h"
#include "libs/sparse.h"

struct MrnaSimulator
{

        // main fields and methods
        ModelParameters* mf_ModelParameters ;
        MrnaSimulator ( ModelParameters* modelParameters , Int seed = 1 ) ;
        void prepareForSteps ( CellState* cellState ) ;
        Doub doStep ( CellState* cellState , Doub targetTime ) ;
        CellState* sampleFromOnlyNativeSteadyState ( Doub toReachSteadyStateDuration ) ;

        // other fields and methods
        Ran ran; // random generator
        Int mm, nn;
        VecDoub a;  // reaction rates
        MatDoub instate, outstate; // reactants matrix, change state matrix
        NRvector<NRsparseCol> outchg, depend; // change state sparse matrix, reaction dependancy sparse matrix
        VecInt pr;  // priority list for reaction
        Doub t; // time
        Doub asum; // sum of all reactions rates
        typedef Doub(MrnaSimulator::*rateptr)(VecDoub &s);
        rateptr *dispatch;
        void describereactions();
        ~MrnaSimulator() {delete [] dispatch;}


	// rate functions
placeholder_rate_functions};




