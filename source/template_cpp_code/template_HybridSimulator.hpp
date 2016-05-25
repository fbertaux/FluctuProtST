/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/

#include "MrnaSimulator.hpp"
#include "HybridRhs.hpp"

# include "libs/odeint.h"

struct HybridSimulator
{

    // constructor and key fields
    HybridSimulator ( ModelParameters* modelParameters , Int randomSeed = 1 ) ;
	ModelParameters* mf_ModelParameters ;
    MrnaSimulator* mf_MrnaSimulator ;

    // key methods
    void simulate ( CellState* cellState , Doub duration ) ;
    void computeGeneMrnaTrajs ( CellState* cellState , Doub duration ) ;


    // mrna/gene trajs storage fields
    const Int MaxEvents ;
    MatDoub FutureMrnaTable;
    MatDoub FutureGeneTable;
    VecDoub FutureGeneMrnaEventsTable;
    Int EventObtained;

    // ode integration fields
    HybridRhs* mf_HybridRhs ;
    Odeint<StepperDopr5<HybridRhs> >* mf_HybridOdeInt ;
    Output* CellOutput ;
    const Doub AbsTolNumErr ;
    const Doub RelTolNumErr ;

};

