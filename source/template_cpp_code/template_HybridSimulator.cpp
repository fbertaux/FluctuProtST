/*
__ FluctuProtST, Version 1.2
__ Francois Bertaux, Inria Paris-Rocquencourt
__ francois.bertaux@inria.fr
__ March 2015
*/

#include "HybridSimulator.hpp"

HybridSimulator::HybridSimulator ( ModelParameters* modelParameters , Int randomSeed ) :
    mf_ModelParameters (modelParameters) , mf_MrnaSimulator ( new MrnaSimulator (mf_ModelParameters,randomSeed) ) ,
    MaxEvents (1e6) , AbsTolNumErr (1e-6) , RelTolNumErr (1e-6)
{
    mf_HybridRhs = new HybridRhs (mf_ModelParameters) ;

    // prepare fields for storing mrna/gene trajs
    FutureGeneTable = MatDoub (mf_ModelParameters->mf_NumGenes,MaxEvents,0.);
    FutureMrnaTable = MatDoub (mf_ModelParameters->mf_NumGenes,MaxEvents,0.);
    FutureGeneMrnaEventsTable = VecDoub (MaxEvents,0.);

    // prepare ode integrator
    CellOutput = new Output () ;
    mf_HybridOdeInt = new Odeint<StepperDopr5<HybridRhs> > ( mf_ModelParameters->mf_NumAllProteinSpecies , 0. , 0. , AbsTolNumErr , RelTolNumErr ,
                                                        0.1 , 0. , *CellOutput , *mf_HybridRhs ) ;
}


void
HybridSimulator::simulate (CellState *cellState, Doub duration)
{
    computeGeneMrnaTrajs ( cellState , duration ) ;
    mf_HybridRhs->setMrnaTables (&FutureMrnaTable,&FutureGeneMrnaEventsTable,EventObtained) ;
    mf_HybridOdeInt->integrate ( cellState->mf_AllProts , 0 , duration ) ;
}


void
HybridSimulator::computeGeneMrnaTrajs (CellState *cellState, Doub duration)
{
    mf_MrnaSimulator->prepareForSteps ( cellState ) ;
    Doub dt,nt; Doub t = 0;
    EventObtained = 1;
    // store init state
    for (int i=0;i<mf_ModelParameters->mf_NumGenes;i++)
    {
        FutureGeneTable[i][0] = cellState->mf_GeneMrnas [3*i] ;
        FutureMrnaTable[i][0] = cellState->mf_GeneMrnas [3*i+2] ;
        FutureGeneMrnaEventsTable[0] = 0. ;
    }
    // simulate
    while (t<duration)
    {
        nt = mf_MrnaSimulator->doStep ( cellState , duration ) ;
        dt = nt - t ;
        if (dt==0) { cout << "dt = 0.." << endl; exit(3);}
        t = mf_MrnaSimulator->t;
        // store
        for (int i=0;i<mf_ModelParameters->mf_NumGenes;i++)
        {
            FutureGeneTable[i][EventObtained] = cellState->mf_GeneMrnas [3*i];
            FutureMrnaTable[i][EventObtained] = cellState->mf_GeneMrnas [3*i+2];
            FutureGeneMrnaEventsTable[EventObtained] = t;
        }
        EventObtained++;
        if (EventObtained>MaxEvents) { cout << "Mrna Table To short. Change Max Events ?" << endl; exit(42); }
    }
}



